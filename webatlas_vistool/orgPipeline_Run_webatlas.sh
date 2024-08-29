##########################
# Pipeline to run webatlas
# Organized by zhuang.liu@ki.se, yongjian.chen@ki.se
# Adapted from https://haniffalab.com/webatlas-pipeline/index.html 

# 1. install the webatlas
wget https://github.com/haniffalab/webatlas-pipeline/archive/refs/tags/v0.5.2.tar.gz
tar -xzvf v0.5.2.tar.gz
cd webatlas-pipeline-0.5.2

# 2. install the nextflow
wget -qO- https://get.nextflow.io | bash
chmod +x nextflow
sudo mv nextflow /usr/local/bin

# 3. install the docker
https://docs.docker.com/desktop/install/mac-install/ 

# 4. select the docker to run the single modality using main.nf function
nextflow run main.nf -params-file /Users/zhuliu/Downloads/Demo/P20063_102_tmp.yaml -entry Full_pipeline -profile docker
nextflow run main.nf -params-file /Users/zhuliu/Downloads/Demo/P20063_102_tmp.yaml -entry Full_pipeline -profile singularity

# 5. using the homebrew (MacOS) to install the npm (Node.js)
brew install node

# 6. run the app with local server
# remember to set the current working path the same as where .json is
cd /Users/zhuliu/Downloads/Demo/temp_Prewebatlas/0.5.2
npx http-server /Users/zhuliu/Downloads/Demo/P20063_102/0.5.1 --port 3000 --cors

# 7. open your Chrome browser 
https://webatlas.cog.sanger.ac.uk/latest/index.html?config=http://localhost:3000/woundhealing-Human-skin-wound-healing-config.json

# 8. Deploy the app the server
......


############################################
# pay attention to the yaml and input files
P20063_102 # put the all the data in one folder
├── 210204_P20063_V10B01-24_LP1_HE_AM_B1-Spot000002.tif # .tif image file
├── analysis # folder from spaceranger output
├── filtered_feature_bc_matrix.h5 # h5 file from spaceranger output
└──spatial # folder from spaceranger output


# Example of P20063_102_tmp.yaml file (10X Visum data)
outdir: /Users/zhuliu/Downloads/Demo/temp_Prewebatlas/

args:
  spaceranger:
    save_h5ad: True

projects:
  - project: woundhealing
    datasets:
      - dataset: Human-skin-wound-healing
        title: "Visium Donor2_Wound1"
        data:
          - data_type: spaceranger
            data_path: /Users/zhuliu/Downloads/Demo/WebAtlas_tools_pre/P20063_102/
          - data_type: raw_image
            data_path: /Users/zhuliu/Downloads/Demo/WebAtlas_tools_pre/P20063_102/210204_P20063_V10B01-24_LP1_HE_AM_B1-Spot000002.tif
          - data_type: label_image_data
            data_path: /Users/zhuliu/Downloads/Demo/WebAtlas_tools_pre/P20063_102/
            file_type: visium
            ref_img: /Users/zhuliu/Downloads/Demo/WebAtlas_tools_pre/P20063_102/210204_P20063_V10B01-24_LP1_HE_AM_B1-Spot000002.tif

vitessce_options:
  spatial:
    xy: "obsm/spatial"
  matrix: "X"
layout: "advanced"



#######################
# Running in the uppmax
cd /crex/proj/snic2021-23-156/publicdata/wound_webatlas/webatlas-pipeline-0.5.2
 
module load bioinfo-tools Nextflow/22.10.1
 
interactive -A naiss2023-22-935 -t 8:00:00 -N 2

nextflow run main.nf -params-file /crex/proj/snic2021-23-156/publicdata/wound_webatlas/input/P20063_102_donor2_wound1_morefun/P20063_102_donor2_wound1_morefun.yaml -entry Full_pipeline -profile singularity
