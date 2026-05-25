# QUALYPSO: Partitioning Uncertainty Components of an Incomplete Ensemble of Climate Projections

## Description
This R package aims at partitioning the total variability in an ensemble of climate projections (or an ensemble of projections of variables impacted by climate change) given that different sources of uncertainty are present. It provides estimates of climate change responses of all simulation chains and of all uncertainty variables. It additionally propagates uncertainty due to missing information in the estimates.

## Getting Started
The main function of the package is simply called "QUALYPSO". At least two arguments must be provided:
- Y is an ensemble of climate projections, i.e. a matrix nS x nY where nS projections are provided for nY years or future time steps;
- scenAvail is a data.frame which provides nEff characteristics for each projection. These characteristics typically indicate the models (GCM, RCM) and emission scenarios (RCP, SSP) that have been used to produce the climate projections. It can also be a spatial entity, or a factor related to an impact model.