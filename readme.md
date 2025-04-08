# An Optimization Framework for Wide-Field Small Aperture Telescope Arrays Used in Sky Surveys

This repository contains the computational framework designed for optimizing wide-field small aperture telescope arrays utilized in sky surveys. Below is a brief description of the key files included in this project:

## File Descriptions

- **violence_search_by_ZEMAX_ori.m**: This is the interactive brute-force search script within the framework. It interacts with ZEMAX to obtain PSF (Point Spread Function) simulated imaging signal-to-noise ratios and calculates the cost of the array configuration.

- **information_collection.m**: This script is designed for collecting training data for the BP neural network based on interactive searches. It gathers data within a certain range corresponding to the telescope and other configurations, which will be used for neural network fitting of the signal-to-noise ratio.

- **train_net_for_each_tel.m**: This script is used for Neural network training for training the corresponding neural network for signal-to-noise ratio calculation.

- **violence_search_by_BP.m**: This is the interactive brute-force search script based on the BP neural network, corresponding to `violence_search_by_ZEMAX_ori.m`.

- **telescopetemplate**: This directory contains the telescope prototype file (.zmx). You can add telescope design files and modify the contents of `violence_search_by_ZEMAX_ori.m` and `violence_search_by_BP.m` accordingly.

- **workingfolder**: This is the working directory where the corresponding working files from the `telescopetemplate` will be copied.

- **00_result**: This folder saves the result files (e.g., figures).
  - The **net_fit_data** folder contains the data files used for training the neural network, as well as the results of the network training.

## Usage

To use this framework, ensure that you have the necessary dependencies installed, including ZEMAX for simulations. Follow the instructions in each script for specific usage details.
