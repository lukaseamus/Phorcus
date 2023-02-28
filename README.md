# Predation governs the eulittoral distribution of a ubiquitous Mediterranean gastropod
This repository contains data and annotated R code accompanying article 10.1007/s10750-023-05143-4 in *Hydrobiologia*. The repository containes the four data and code folders **Introduction**, **Density**, **Predation** and **Physiology** alongside all final **Figures** in vector format. Below is a description of each file within the data and code folders as well as the input and output of each R script.

**Introduction**
1. `Steromphala.csv`: Data on the lower and upper limits of the distributions of Steromphala cineraria and S. umbilicalis.
    - **species** = categorical variable with levels *Laminaria digitata* (d), *Laminaria hyperborea* (h) and *Laminaria ochroleuca* (o)
    - **bag** = random factor (categorical variable) with plants (P3...24) and mesh bags (B1...9) as levels
    - **age** = detrital age given in days
    - **R** = respiration rate given in µmol oxygen per gram of buoyant mass per hour (buoyant mass is practically identical to wet mass)
    - **NPP** = net photosynthesis rate given in µmol oxygen per gram of buoyant mass per hour
    - **GPP** = gross photosynthesis rate given in µmol oxygen per gram of buoyant mass per hour (NPP + R)
    - **d:w** = dry to wet mass ratio
2. `Steromphala.R`: Code to analyse and visualise carbon assimilation.
    - **Input** = `Assimilation.csv`, Figure 5b from `Decomposition.R`
    - **Output** = Figure 5, Figure S8, carbon assimilation functions and results

**Density**
1. `Export.csv`: Carbon export data.
    - **month** = month and year given as MMM-YY
    - **season** = categorical variable with levels Spring, Summer, Autumn and Winter
    - **time** = numerical expression of months
    - **species** = categorical variable with levels *Laminaria digitata* (d), *Laminaria hyperborea* (h) and *Laminaria ochroleuca* (o)
    - **dw.export** = biomass export given in grams of dry mass per plant per day
    - **fw.export** = biomass export given in grams of wet mass per plant per day, converted from dry mass with plant-specific dry to wet mass ratios
    - **fw.export.avg** = biomass export given in grams of wet mass per plant per day, converted from dry mass with species- and month-specific dry to wet mass ratios
    - **C.export** = carbon export given in grams per plant per day, converted from dry mass with species- and month-specific carbon content (%)
2. `Carbon.csv`: Lamina carbon content data.
    - **month** = month and year given as MMM-YY
    - **season** = categorical variable with levels Spring, Summer, Autumn and Winter
    - **species** = categorical variable with levels *Laminaria digitata* (d), *Laminaria hyperborea* (h) and *Laminaria ochroleuca* (o)
    - **carbon** = carbon content (%)
3. `Mass.csv`: Sporophyte mass data.
    - **month** = month and year given as MMM-YY
    - **season** = categorical variable with levels Spring, Summer, Autumn and Winter
    - **species** = categorical variable with levels *Laminaria digitata* (d), *Laminaria hyperborea* (h) and *Laminaria ochroleuca* (o)
    - **mass** = whole plant wet mass given in grams
4. `DW.csv`: Dry to wet mass ratio data.
    - **month** = month and year given as MMM-YY
    - **season** = categorical variable with levels Spring, Summer, Autumn and Winter
    - **species** = categorical variable with levels *Laminaria digitata* (d), *Laminaria hyperborea* (h) and *Laminaria ochroleuca* (o)
    - **d.w** = dry to wet mass ratio
5. `Density.csv`: Sporophyte density data.
    - **month** = month and year given as MMM-YY
    - **season** = categorical variable with levels Spring, Summer, Autumn and Winter
    - **species** = categorical variable with levels *Laminaria digitata* (d), *Laminaria hyperborea* (h) and *Laminaria ochroleuca* (o)
    - **density** = number of plants per square metre. Note that *Laminaria hyperborea* and *Laminaria ochroleuca* share each quadrat while *Laminaria digitata* occurs spatially separated. 
6. `Export.R`: Code to analyse and visualise carbon export.
    - **Input** = `Export.csv`, `Carbon.csv`, `Mass.csv`, `DW.csv`, `Density.csv`
    - **Output** = Figure S2, `Constants.csv`, carbon export results

**Predation**
1. `Decomposition.csv`: Biomass decomposition data.
    - **species** = categorical variable with levels *Laminaria digitata* (d), *Laminaria hyperborea* (h) and *Laminaria ochroleuca* (o)
    - **site** = categorical variable with levels West Hoe (WH), Drake's Island (DI) and Jennycliff (JC)
    - **substratum** = categorical variable with levels Forest and Sediment
    - **mesh** = mesh diameter given in centimetres
    - **g.loss** = absolute biomass loss given in grams per day
    - **perc.loss** = relative biomass loss given in percentage of initial mass per day
2. `Biochemical.csv`: Elemental stoichiometry and phenols in relation to decomposition.
    - **species** = categorical variable with levels *Laminaria digitata* (d), *Laminaria hyperborea* (h) and *Laminaria ochroleuca* (o)
    - **bag** = random factor (categorical variable) with mesh bags (B1...9) as levels
    - **age** = detrital age given in days
    - **g.loss** = absolute biomass loss given in grams per day
    - **perc.loss** = relative biomass loss given in percentage of initial mass per day
    - **phenols** = final soluble polyphenolic content (%)
    - **N** = final nitrogen content (%)
    - **C** = final carbon content (%)
    - **CN** = final carbon to nitrogen ratio
3. `Grazing.csv`: Image analysis data of tissue damage on final retrieval.
    - **species** = categorical variable with levels *Laminaria digitata* (d), *Laminaria hyperborea* (h) and *Laminaria ochroleuca* (o)
    - **bag** = random factor (categorical variable) with mesh bags (B1...9) as levels
    - **excavation** = surface area of excavation scars relative to total tissue area (%)
    - **perforation** = surface area of holes relative to total tissue area plus holes (%)
4. `Decomposition.R`: Code to analyse and visualise carbon export.
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
