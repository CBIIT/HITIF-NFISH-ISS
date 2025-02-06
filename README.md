# Project Name

This repository contains [Brief Project Description]. Follow the instructions below to set up and run the project.

## **Installation Guide**

### **Step 1: Install Anaconda**
Anaconda is a Python distribution that includes package management and environment handling.

1. **Download Anaconda:**
   - Go to [Anaconda Official Website](https://www.anaconda.com/products/distribution)
   - Download the appropriate installer for your OS (Windows, macOS, or Linux)
   - Follow the installation instructions provided on the website

2. **Verify Installation (Optional)**
   After installation, open a terminal (Command Prompt, PowerShell, or Terminal) and run:

   ```bash
   conda --version
    ```
   



### **Step 2: Clone This Repository**

Once Anaconda is installed, clone this repository to your local machine.

#### **Using Git (Recommended)**
If you have Git installed, open a terminal and run:

```bash
git clone https://github.com/CBIIT/HITIF-NFISH-ISS.git
cd HITIF-NFISH-ISS
```


### **Step 3: Create a new conda environment**

Create a new conda environment using the requirements.txt file:

1. Create a new Conda environment with a specified Python version

Open a terminal and run:
```
conda create --name hitif-nfish-ops python=3.10
```

2. Activate the new environment
```
conda activate hitif-nfish-ops
```
3. Install dependencies from the requirements file
```
conda install --file requirements_nfish.txt
```




