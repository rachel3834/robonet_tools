# robonet_tools
Suite of useful tools developed in the course of operating the ROME/REA Project

## Installation

This package consists of a number of tools designed to be run manually as needed, usually from an interactive commandline. 
To install it, simply clone this repository to a convenient location on your local machine

## Configuration

For tools with many configurable parameters, this package uses JSON format configuration files.  Template JSON configuration files can be found in the configs/ sub-directory.  To configure your local installation it is recommended that you copy the template configuration files to a different location that the cloned repository on your local machine, and then edit the configuration files to reflect your desired configuration.  

The local, customized copy of the configuration file is then generally called as a commandline argument for scripts that require the configuration.  

For example, to configure and run lightcurve_download_aws.json:
* Copy the robonet_tools/configs/lightcurve_download_aws.json file to a different location and update this file to reflect your AWS user ID and local paths.  This is the version of this file that you should call when you run the download script.  

* To run the script, type [Unix, MacOS]:
> python download_lightcurves_aws.py /path/to/local/lightcurve_download_aws.json

 or in Windows:

> C:\path\to\python\installation\python.exe "C:/path/to/cloned/repository/robonet_tools/data_harvest/download_lightcurves_aws.py /path/to/local/lightcurve_download_aws.json"
