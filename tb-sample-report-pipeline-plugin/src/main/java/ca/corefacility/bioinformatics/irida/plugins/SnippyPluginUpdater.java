package ca.corefacility.bioinformatics.irida.plugins;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Path;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.Scanner;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import com.google.common.base.Joiner;
import com.google.common.base.Splitter;
import com.google.common.collect.ImmutableMap;
import com.google.common.collect.Lists;
import com.google.common.collect.Sets;

import ca.corefacility.bioinformatics.irida.exceptions.IridaWorkflowNotFoundException;
import ca.corefacility.bioinformatics.irida.exceptions.PostProcessingException;
import ca.corefacility.bioinformatics.irida.model.sample.MetadataTemplateField;
import ca.corefacility.bioinformatics.irida.model.sample.Sample;
import ca.corefacility.bioinformatics.irida.model.sample.metadata.MetadataEntry;
import ca.corefacility.bioinformatics.irida.model.sample.metadata.PipelineProvidedMetadataEntry;
import ca.corefacility.bioinformatics.irida.model.workflow.IridaWorkflow;
import ca.corefacility.bioinformatics.irida.model.workflow.analysis.AnalysisOutputFile;
import ca.corefacility.bioinformatics.irida.model.workflow.analysis.type.AnalysisType;
import ca.corefacility.bioinformatics.irida.model.workflow.submission.AnalysisSubmission;
import ca.corefacility.bioinformatics.irida.pipeline.results.updater.AnalysisSampleUpdater;
import ca.corefacility.bioinformatics.irida.service.sample.MetadataTemplateService;
import ca.corefacility.bioinformatics.irida.service.sample.SampleService;
import ca.corefacility.bioinformatics.irida.service.workflow.IridaWorkflowsService;

import com.fasterxml.jackson.databind.ObjectMapper;
import com.fasterxml.jackson.core.type.TypeReference;

import java.util.List;

import com.fasterxml.jackson.annotation.JsonProperty;
import com.fasterxml.jackson.annotation.JsonIgnoreProperties;

class QC {
    public float pct_reads_mapped;
    public float num_reads_mapped;

}

class Lineage {
    public String lin;
    public String family;
    public String spoligotype;
    public String rd;
    public float frac;
}

@JsonIgnoreProperties(ignoreUnknown = true)
class Variant {
    public String sample;
    public String gene_name;
    public String chr;
    public int genome_pos;
    public String type;
    public String change;
    public float freq;
    public String nucleotide_change;
    public String locus_tag;
    public String gene;
    public String drug;
    public String confidence;
}

class DBVersion {
    public String name;
    public String commit;
    public String author;
    public String date;

    @JsonProperty("Author")
    public String getAuthor() {
        return author;
    }

    @JsonProperty("Author")
    public void setAuthor(String author) {
        this.author = author;
    }

    @JsonProperty("Date")
    public String getDate() {
        return date;
    }

    @JsonProperty("Date")
    public void setDate(String date) {
        this.date = date;
    }
}

class Pipeline {
    public String mapper;
    public String variant_caller;
}

@JsonIgnoreProperties(ignoreUnknown = true)
public class TBProfilerReport {
    public QC qc;
    // region_coverage key not processed
    public List<Lineage> lineage;
    public String main_lin;
    public String sublin;
    public List<Variant> dr_variants;
    public List<Variant> other_variants;
    public String XDR;
    public String MDR;
    public String drtype;
    public DBVersion db_version;
    public String id;
    public String tbprofiler_version;
    public Pipeline pipeline;
}

class MetadataValue {
	public String header;
	public String value;

	MetadataValue(String header, String value) {
		this.header = header;
		this.value = value;
	}
}

/**
 * This implements a class used to perform post-processing on the analysis
 * pipeline results to extract information to write into the IRIDA metadata
 * tables. Please see
 * <https://github.com/phac-nml/irida/blob/development/src/main/java/ca/corefacility/bioinformatics/irida/pipeline/results/AnalysisSampleUpdater.java>
 * or the README.md file in this project for more details.
 */
public class SnippyPluginUpdater implements AnalysisSampleUpdater {
	private static final Logger logger = LoggerFactory.getLogger(SnippyPluginUpdater.class);

	private static final String TBPROFILER_FILE = "tb_profiler_json.json"; 

	private final MetadataTemplateService metadataTemplateService;
	private final SampleService sampleService;
	private final IridaWorkflowsService iridaWorkflowsService;

	// @formatter:off
	private Map<String, String> TBPROFILER_FIELDS = ImmutableMap.<String,MetadataValue>builder()
		.put("drtype", new MetadataValue("Drug Resistance Type", ""))
		.put("lineage", new MetadataValue("Lineage", ""))
		.put("family", new MetadataValue("Family", ""))
		.put("spoligotype", new MetadataValue("Spoligotype", ""))
		.put("lineage_agreement", new MetadataValue("% Lineage Agreement", ""))
		.put("isoniazid", new MetadataValue("Isoniazid", ""))
		.put("rifampicin", new MetadataValue("Rifampicin", ""))
		.put("ethambutol", new MetadataValue("Ethambutol", ""))
		.put("streptomycin", new MetadataValue("Streptomycin", ""))
		.put("other_resistance", new MetadataValue("Other resistance", ""))
		.put("isoniazid_variants", new MetadataValue("Isoniazid Res Variants", ""))
		.put("rifampicin_variants", new MetadataValue("Rifampicin Res Variants", ""))
		.put("ethambutol_variants", new MetadataValue("Ethambutol Res Variants", ""))
		.put("streptomycin_variants", new MetadataValue("Streptomycin Res Variants", ""))
		.put("other_resistance_variants", new MetadataValue("Other resistance Res Variants", ""))
		.put("tbprofiler_version", new MetadataValue("TbProfiler Version", ""))
		.build();
			// @formatter:on


	/**
	 * Builds a new {@link SnippyPluginUpdater} with the given services.
	 * 
	 * @param metadataTemplateService The metadata template service.
	 * @param sampleService           The sample service.
	 * @param iridaWorkflowsService   The irida workflows service.
	 */
	public SnippyPluginUpdater(MetadataTemplateService metadataTemplateService, SampleService sampleService,
			IridaWorkflowsService iridaWorkflowsService) {
		this.metadataTemplateService = metadataTemplateService;
		this.sampleService = sampleService;
		this.iridaWorkflowsService = iridaWorkflowsService;
	}

	/**
	 * Code to perform the actual update of the {@link Sample}s passed in the
	 * collection.
	 * 
	 * @param samples  A collection of {@link Sample}s that were passed to this
	 *                 pipeline.
	 * @param analysis The {@link AnalysisSubmission} object corresponding to this
	 *                 analysis pipeline.
	 */
	@Override
	public void update(Collection<Sample> samples, AnalysisSubmission analysis) throws PostProcessingException {
		if (samples == null) {
			throw new IllegalArgumentException("samples is null");
		} else if (analysis == null) {
			throw new IllegalArgumentException("analysis is null");
		} else if (samples.size() != 1) {
			// In this particular pipeline, only one sample should be run at a time so I
			// verify that the collection of samples I get has only 1 sample
			throw new IllegalArgumentException(
				"Expected one sample; got '" + "samples size=" + samples.size() + " is not 1 for analysisSubmission=" + analysis.getId());
		}

		// extract the 1 and only sample (if more than 1, would have thrown an exception
		// above)
		final Sample sample = samples.iterator().next();

		// extracts paths to the analysis result files
		AnalysisOutputFile tbprofilerFile = analysis.getAnalysis().getAnalysisOutputFile(TBPROFILER_FILE);
		Path filePath = tbprofilerFile.getFile();

		try {
			Map<String, MetadataEntry> metadataEntries = new HashMap<>();

			// get information about the workflow (e.g., version and name)
			IridaWorkflow iridaWorkflow = iridaWorkflowsService.getIridaWorkflow(analysis.getWorkflowId());
			String workflowVersion = iridaWorkflow.getWorkflowDescription().getVersion();
			String workflowName = iridaWorkflow.getWorkflowDescription().getName();

			//Read the JSON file from TbProfile output
			@SuppressWarnings("resource")
			String jsonFile = new Scanner(new BufferedReader(new FileReader(filePath.toFile()))).useDelimiter("\\Z").next();

			// map the results into a Map
			ObjectMapper mapper = new ObjectMapper();
			TBProfilerReport tbprofilerResults = mapper.readValue(jsonFile, new TypeReference<Map<String, Object>>() {});

			TBPROFILER_FIELDS.entrySet()
			if (tbprofilerResults.size() > 0) {
				//Map<String, Object> result = tbprofilerResults.get(0);
				tbprofilerResults.dr_variants.forEach(variant -> {
					if (TBPROFILER_FIELDS.containsKey(variant.drug)) {
						String drug_resistant_variants = ()
					}
				});
				

				//loop through each of the requested fields and append workflow version and save the entries
				TBPROFILER_FIELDS.entrySet().forEach(e -> {
					if (tbprofilerResults.containsKey(e.getKey())) {
						Object valueObject = tbprofilerResults.get(e.getKey());
						String value = (valueObject != null ? valueObject.toString() : "");
						PipelineProvidedMetadataEntry metadataEntry =
								new PipelineProvidedMetadataEntry(value, "text", analysis);
						metadataEntries.put(e.getValue() + " (v"+workflowVersion+")", metadataEntry);
					}
				});

				// convert string map into metadata fields
				Map<MetadataTemplateField, MetadataEntry> metadataMap = metadataTemplateService.getMetadataMap(metadataEntries);

				//save metadata back to sample
				samples.forEach(s -> {
					s.mergeMetadata(metadataMap);
					sampleService.updateFields(s.getId(), ImmutableMap.of("metadata", s.getMetadata()));
				});

			} else {
				throw new PostProcessingException("TbProfiler results for file are not correctly formatted");
			}
			
		} catch (IOException e) {
			throw new PostProcessingException("Error parsing JSON from TbProfiler (Snippy-Tb-Sample-Report) results", e);
		} catch (IridaWorkflowNotFoundException e) {
			throw new PostProcessingException("Could not find workflow for id=" + analysis.getWorkflowId(), e);
		}
	}
	
	/**
	 * The {@link AnalysisType} this {@link AnalysisSampleUpdater} corresponds to.
	 * 
	 * @return The {@link AnalysisType} this {@link AnalysisSampleUpdater}
	 *         corresponds to.
	 */
	@Override
	public AnalysisType getAnalysisType() {
		return SnippyPipelinePlugin.VARIANT_CALLING;
	}
}
