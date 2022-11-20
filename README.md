![Title](Banner/Solar_Analysis_Toolkit.png)

Python library combining multiple solar analysis methods into one convenient package.

Features a variety of methodologies commonly used to analyse features, events, and statistics of solar data.

As of now, only Windows is currently supported. Pull requests are being made daily to support use on Linux distributions.

# **Installation**
TODO

For now, users will have to download the repository and solve for the dependencies on their own. Works best with Python Anaconda, and it is recommended to start off in a clean enviornment 


# **Usage**
Several programs are built directly within the main directory of Solar_Analysis_Toolkit. To use these programs, just simply call the program name and input data when prompted. However, some programs are dependent on others (e.g. Run_Solar_MFDFA_Analysis.py depends on already running Run_Solar_MFDFA.py, which also requires a time series data set to be obtained, which can be created using Create_AIA_Time_Series.py). For now, one can use the toolkit to identify events that happen on the sun which are captured by the aforementioned SDO spacecraft. Currently, the only telescope of interest that is supported is the AIA/HMI instruments onboard the NASA SDO spacecraft. Future goals look to improve and implement more analysis methods regarding different spacecrafts and missions such as the Parker Solar Probe (PSP) or HINODE.

