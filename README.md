# Monitoring Water Hyacinth using EO in Winam Gulf, Lake Victoria (2018-2021)





<img width="974" alt="Screenshot 2025-03-24 at 10 55 52" src="https://github.com/user-attachments/assets/6ea09a59-07c7-42d9-8314-6ffb289798dc" />



# Background
Water hyacinth is considered one of the world's worst invasive aquatic weeds. It grows in thick, free-floating mats that can double in extent in 5-15 days (see figure 1). Native to the Amazon, the plant was introduced in Lake Victoria in 1989 as an ornamental pond plant, causing periodic waves of infestation bearing socioeconomic and environmental issues. Previous studies monitored water hyacinth invasions in Lake Victoria to inform management practices until 2017. This study proposes a Google Earth Engine EO framework to monitor the weed in Winam Gulf (Kenya), the area of the Lake most historically affected by water hyacinth invasions, over the period 2018-2021 using both active and passive remote sensing.



# Aims
- Develop Google Earth Engine scripts to classify and measure the extent of water hyacinth in Winam Gulf from Sentinel-2 and Sentinel-1 imagery over the period 2018-2021.
- Identify trends from Sentinel-2 and Sentinel-1 2018-2021 water hyacinth extent estimates and compare them to weather data from Kisumu, Kenya.
- Create Sentinel-1 yearly heatmaps to identify areas in the Gulf that have been most prone to water hyacinth invasion in the period 2018-2021.



# EO methods
Google Earth Engine workflows to classify water haycinth were developed for both Sentinel-2 and Sentinel-2 data. The figure below shows the sampling points used in this study for Random Forest classification using Sentinel-2 data. Due to Sentinel-2 L2A data unavailability in the first 11 months of 2018, Sentinel-2 L1C images were manually processed and thresholded to cover this period. A Sentinel-1 script was also developed to create and threshold biweekly radar composites, as well as for creating yearly annual heatmaps of water hyacinth invasion. The classification scripts were used to obtain a timeseries of the weed's extent over 2018-2021. 

<img width="443" alt="Screenshot 2025-03-24 at 10 58 28" src="https://github.com/user-attachments/assets/de65db85-749f-4f0b-bdd6-a54d3a646a4c" />


# To see the study's findings
visit the study's website available at: https://www.geos.ed.ac.uk/~mscgis/21-22/s1737656/



