# CGA-BNI
CGA-BNI is a promising and scalable tool to predict both the structure and the dynamics of a gene regulatory network when a highest accuracy is needed at the cost of sacrificing the execution time.

# Guide for executing the tool
## Setup the tool

- Install Java SE Development Kit, recommended version is 8 (from https://www.oracle.com/java/technologies/javase/javase-jdk8-downloads.html). We are using the java version "1.8.0_171"
- Install NetBeans IDE, recommeded version is above 7.3 (from https://www.oracle.com/technetwork/java/javase/downloads/jdk-netbeans-jsp-3413139-esa.html)
- Download this project into a folder, for example, "**CGAProj**"
- Extract all **.zip files** in the **/data/Inference** folder
- Open the project by NetBeans IDE, and open the file "**InferBN.java**":
  - The **main** function, and also inference funtions is in this Java file
  - Build and run the project with the **default config**:
    - The config has **three arguments**: the path to data, type of data, number of trials
    
      + **path to data**: the path to the data folder containing DREAM3, random networks, and the large-scale E.coli datasets. We could make another dataset by making another dataset folder and chaning the source codes a little similar to the way of the large-scale Ecoli dataset.
      + **type of data**: **"dream"**, **"rbn"**, and **"ecoli"** correspond to DREAM3, random networks, and the large-scale Ecoli dataset, respectively
      + **number of trials**: repeat the experiment a number of times to get avergage results

    - Right-click in the project, choose menu **Properties -> Run**
      ![Project default configuration](/images/cga1.png)
      
      + We could change the arguments as in the above figure

    - Some examples of arguments:
      + E:\\CGA_BNI\\data\\ dream 1
      + E:\\CGA_BNI\\data\\ rbn 50
      + E:\\CGA_BNI\\data\\ ecoli 20
      
- Results are stored in the project folder and named as **run1**, **run2**, ...

## Contacts 

Contact us: trinhhungcuong@tdtu.edu.vn

Thanks for your interests!


