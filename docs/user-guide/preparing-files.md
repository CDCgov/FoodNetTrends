## Preparing Files

Data preparation is crucial for accurate modeling. The FoodNet enhanced model requires the following inputs:

- **MMWR Data File** (`mmwrFile`): The MMWR (Morbidity and Mortality Weekly Report) data file in SAS format (`.sas7bdat`).
- **Census Data Files**: Two census data files in SAS format:
  - **Bacterial Census File** (`censusFile_B`)
  - **Parasitic Census File** (`censusFile_P`)
- **Parameters**:
  - **Travel Status** (`travel`): A list of travel statuses to include (e.g., `"NO,UNKNOWN,YES"`).
  - **CIDT Variables** (`cidt`): A list of CIDT (Culture-Independent Diagnostic Tests) variables (e.g., `"CIDT+,CX+,PARASITIC"`).
  - **Project ID** (`projID`): A unique project identifier (e.g., `"20240706"`).

Ensure that the data files are accessible and properly formatted.
