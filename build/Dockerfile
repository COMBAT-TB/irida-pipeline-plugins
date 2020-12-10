FROM ubuntu:20.04 as builder
MAINTAINER Peter van Heusden <pvh@sanbi.ac.za>

ENV IRIDA_TAG=20.09.3
RUN apt update && apt install -y default-jdk-headless maven git
RUN git clone --depth 1 -b $IRIDA_TAG https://github.com/phac-nml/irida.git && cd irida/lib && bash install-libs.sh && cd ..
RUN cd irida && mvn clean install -DskipTests

FROM alpine:latest
RUN apk update && apk add openjdk11-jre-headless maven
COPY --from=builder /root/.m2 /root/
RUN mkdir /root/.m2 && mv /root/repository /root/.m2/
