# Predation governs the eulittoral distribution of a ubiquitous Mediterranean gastropod
This repository contains data and annotated R code accompanying article 10.1007/s10750-023-05143-4 in *Hydrobiologia*. The repository containes the four data and code folders **Introduction**, **Density**, **Predation** and **Physiology** alongside all final **Figures** in vector format. Below is a description of each column within the data and the input and output of each R script.

**Introduction**
1. `Steromphala.csv`: Data on the lower and upper limits of the distributions of *Steromphala cineraria* and *S. umbilicalis*.
    - **reference** = author name and year of referenced data
    - **DOI** = digital object identifier or uniform resource locator (URL) of referenced data
    - **species** = categorical variable with levels *Steromphala cineraria* and *Steromphala umbilicalis*
    - **lower** = lower distributional limit in relation to lowest astronomical tide given in metres
    - **upper** = upper distributional limit in relation to lowest astronomical tide given in metres
2. `Steromphala.R`: Code to calculate introductory descriptive statistics.
    - **Input** = `Steromphala.csv`
    - **Output** = descriptive statistics in paragraph two of the introduction

**Density**
1. `density.csv`: Density data for *Phorcus turbinatus*, *Stramonita haemastoma* and *Thalassoma pavo*.
    - **site** = categorical variable with levels Xwejni, Dwejra and Ras 
    - **date** = date given as DD.MM.YY
    - **species** = categorical variable with levels *Phorcus turbinatus*, *Stramonita haemastoma* and *Thalassoma pavo*
    - **original** = count of individuals in defined area (40-cm quadrat for *Phorcus turbinatus* and *Stramonita haemastoma*, 1×10-m transect for *Thalassoma pavo*)
    - **adjusted** = density given in individuals per square metre
2. `density.R`: Code to analyse and visualise *Phorcus turbinatus* density and distribution data.
    - **Input** = `density.csv`, `distribution.csv` from **Predation**
    - **Output** = Figure 1

**Predation**
1. `distribution.csv`: Distribution data for all studied species.
    - **site** = categorical variable with levels Xwejni, Dwejra and Ras
    - **date** = date given as DD.MM.YY
    - **time** = time given as HH:MM
    - **species** = categorical variable with levels *Phorcus turbinatus*, *Stramonita haemastoma*, *Thalassoma pavo*, *Hermodice carunculata* and *Hexaplex trunculus*
    - **position** = position in metres in relation to sea level at the time of observation
    - **tide** = tidal level in metres above lowest astronomical tide at the time of observation provided by the Lampedusa tide station 
    - **lat.position** = position in relation to lowest astronomical tide given in metres
2. `predation.csv`: *In situ* predation data.
    - **site** = categorical variable with levels *Laminaria digitata* (d), *Laminaria hyperborea* (h) and *Laminaria ochroleuca* (o)
    - **date** = random factor (categorical variable) with mesh bags (B1...9) as levels
    - **time** = detrital age given in days
    - **position** = absolute biomass loss given in grams per day
    - **tide** = relative biomass loss given in percentage of initial mass per day
    - **lat.position** = final soluble polyphenolic content (%)
    - **length** = final nitrogen content (%)
    - **mass** = final carbon content (%)
    - **predation** = final carbon to nitrogen ratio
3. `muricids.csv`: *In vitro* predation data.
    - **species** = categorical variable with levels *Laminaria digitata* (d), *Laminaria hyperborea* (h) and *Laminaria ochroleuca* (o)
    - **bag** = random factor (categorical variable) with mesh bags (B1...9) as levels
    - **excavation** = surface area of excavation scars relative to total tissue area (%)
    - **perforation** = surface area of holes relative to total tissue area plus holes (%)
4. `predation.R`: Code to analyse and visualise carbon export.
    - **Input** = `Decomposition.csv`, `Biochemical.csv`, `Grazing.csv` 
    - **Output** = Figure 3, Figure S3, decomposition results
5. `muricids.R`: Code to analyse and visualise carbon export.
    - **Input** = `Decomposition.csv`, `Biochemical.csv`, `Grazing.csv` 
    - **Output** = Figure 3, Figure S3, decomposition results
    
**Physiology**
1. `L4.csv`: Physical and chemical data from station L4, compiled from data available at https://www.westernchannelobservatory.org.uk/l4_ctdf/index.php.
    - **Date** = date given as DD.M.YY
    - **Month** = month
    - **Year** = year given as YYYY
    - **Season** = categorical variable with levels Spring, Summer, Autumn and Winter
    - **Temp** = temperature (°C)
    - **Fluor** = fluorescence given in milligrams of chlorophyll *a* per cubic metre
    - **Depth** = depth given in metres
    - **Density** = water density given in kilograms per cubic metre
    - **Salinity** = salinity (‰)
    - **Trans** = transmission (%)
    - **PAR** = photosynthetically active radiation given in µmol photons per square metre per second
    - **Oxygen** = oxygen given in µM
    - **Sound** = sound velocity given in metres per second
2. `Irradiance.R`: Code to analyse the depth-irradiance relationship.
    - **Input** = `L4.csv`
    - **Output** = seasonal and annual exponential depth-irradiance relationships 

Luka Seamus Wright, 28 February 2023
