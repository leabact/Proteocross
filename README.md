# Proteocross
Proteocross - proteomics data treatment automation

This app is made to help you (and your colleagues) handle the results from two expriments of protein identification : cutting and digestion of SDS-PAGE gel bands after migration, and pull-down (capture and enrichment of target.s protein.s). You might want to reduce the Excel use to a minimum : in that case, this app could give you a hand !

### This code will help you if :
- You cut several bands from your protein migration gel and want to know if they share proteins of interest
- You made a pull-down to enrich your protein.s of interest and want to know if they match the ones identified in the cut band.s

It crosses the results obtained from the identification of proteins from the two techniques.

### About the band.s
- You will be asked how many bands were cut and to select as many excel file.s containing the identified proteins
- For each band, you will have to enter the minimal and maximal molecular weight of the protein of interest

### About the pull-down
- You will be asked to enter the lowest ratio that you consider significant for your experiment
- A protein will be considered "of interest" if its ratio is > than the one you entered and its p-value ≤ 0.05

If you have pull-down data, the app will plot the associated volcano, highlighting the proteins shared by the pull down and each band. You will also have a 'Top 15 Ratio' volcano plot and associated data : it highlights the 15 significants proteins having the best ratio.

### Results
Results will be Excel files -if you have pull-down data,associated volcanos plot will be .png- and will be saved in the same folder as your data.

You can run as many analyses as you want or need one after the other, but make sure to have your different data in separate folders, otherwise the new Excel results files will overwrite the previous ones.
