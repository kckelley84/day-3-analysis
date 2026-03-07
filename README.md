# Day-3 Assisted Hatching and Embryo Development
## Overview

Assisted hatching (AH) is sometimes performed on Day-3 embryos to facilitate later blastocyst expansion. However, its effects on embryo development and genetic outcomes remain unclear.

This project evaluates whether Day-3 assisted hatching influences embryo development, morphology, or ploidy status using time-lapse embryo imaging and mixed-effects statistical models.

The analysis examines associations between assisted hatching and:

- AI-derived embryo quality score

- blastocyst morphology (Gardner scores)

- developmental timing

- blastulation probability

- embryo euploidy

## Dataset

Embryos from a retrospective IVF dataset were analyzed.

Study population

- 470 embryos reaching t8

- 351 biopsied embryos

- Mean patient age: 33.1 years

Development outcomes
|Outcome|	Count|
|----|---|
|Biopsied embryos	|351|
|Euploid	|207|
|Aneuploid	|128|
|Mosaic|	10|
|No result|	6|

Two embryo groups were compared:

|Group|	Embryos|
|----|---|
|Day-3 Assisted Hatching|	220|
|No Assisted Hatching	|250|

Mean EQ scores were similar between groups.

## Analytical Approach

Mixed effects regression models were used to evaluate associations between assisted hatching and embryo outcomes.

All models included:

Predictors

- Day-3 assisted hatching

- maternal age

Random effect

- patient ID

This accounts for multiple embryos per patient.

## Results
### AI Embryo Quality Score

No significant association was observed between assisted hatching and AI score.

|Metric|	Value|
|----|---|
|β	|0.168|
|SE	|0.123|
|p-value|	0.17|

### Blastocyst Morphology
Expansion stage

Day-3 assisted hatching was associated with higher blastocyst expansion stage.

|Metric|	Value|
|----|---|
|β	|0.706|
|p-value|	<0.001|

This corresponds to an average increase of ~0.7 expansion stage.

ICM and TE quality

No significant association was observed.

|Outcome|	p-value|
|----|---|
|ICM quality	|0.278|
|TE quality	|0.609|

Developmental Timing

Assisted hatching delayed the time to expanded blastocyst.

|Metric|	Value|
|----|---|
|Mean tEB (AH)|	122.6 hours|
|Mean tEB (No AH)|	117.1 hours|
|Difference|	+5.5 hours|
|p-value	|<0.001|

Blastulation Probability

Assisted hatching was not associated with probability of reaching tEB.

|Group	|Blastulation Rate|
|----|---|
|Day-3 AH	|65.4%|
|No AH	|60.8%|
|p-value| 0.193|

Euploidy

No association between assisted hatching and euploidy.

|Metric	|Value|
|----|---|
|OR|	0.904|
|95% CI|	0.524–1.560|
|p-value	|0.717|

## Key Findings

- Assisted hatching increases blastocyst expansion stage

- Assisted hatching delays time to expanded blastocyst by ~5 hours

- Assisted hatching does not affect embryo quality scores

- Assisted hatching does not affect blastulation probability

- Assisted hatching does not affect euploidy

Together these findings suggest that assisted hatching influences blastocyst morphology but not embryo viability.

## Skills Demonstrated

- mixed effects regression modeling

- hierarchical biomedical datasets

- embryo morphokinetic analysis

- clinical intervention evaluation

- interpretation of biological timing data

## Data Privacy

Patient identifiers and proprietary clinical data have been removed. Example or simulated data may be included for demonstration purposes.
