# Simulation Code
The code used for simulating the ultralight-scalar/black hole interactions, and processing it for data analysis (presentation of GW strains, drifts, etc.)

## Workflow
An example workflow for running simulations over a selection of input parameters is given below. Each simulation (step 3) is independent, so can (and in general, should) be run in parallel.
1. Modify the config files to replicate the conditions you want to simulate. Note that directories are not created on the fly; if you want a `MOA/results` folder, this must be manually created. By default, the `Config.py` file is configured to generate a small example for MOA-2011. The required directories are present by default in the repository.
2. Create the input files for the desired configuration:
```python ConditionMaker.py```
3. Run the simulation taking each of the input fies as input. In our example:  
```python ULBSimulator.py exampleMOA/inputs/0X0.txt```  
```python ULBSimulator.py exampleMOA/inputs/0X1.txt```  
```python ULBSimulator.py exampleMOA/inputs/1X0.txt```  
```python ULBSimulator.py exampleMOA/inputs/1X1.txt```  
4. Run the file reader over each of the output files. In our example:  
```python FileReader.py exampleMOA/outputs/0X0M7.3E+00a9.9E-01f1.0E+12mu1.8E-13.txt```  
```python FileReader.py exampleMOA/outputs/0X1M7.3E+00a9.9E-01f1.0E+12mu3.7E-12.txt```  
```python FileReader.py exampleMOA/outputs/1X0M7.3E+00a9.9E-01f1.0E+20mu1.8E-13.txt```  
```python FileReader.py exampleMOA/outputs/1X1M7.3E+00a9.9E-01f1.0E+20mu3.7E-12.txt```
5. The files in the results folder are now processed and ready for analysis. 
