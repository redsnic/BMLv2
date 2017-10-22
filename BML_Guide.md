# Guide to BML graphical user interface

## Introduction
**BML** is a software to inferr the pathways of cancer evolution from the data of **cross-sectional** studies cerated by Navodit Mistra et al., for more informations see **'BML publication'** in Help menu.
## Input formats
**BML** natively supports only a custom format based on a **boolean matrix** that associates **samples** (on the y-axis) and **genes** (on the x-axis) placing a **1** for each sample that shows a mutation in the given gene, and **0** otherwise. The format is so organized: 

* The first line contains two numbers, the **total number of samples** and the **total number of genes**, spaced by a tab.
* The second line contains the **Hugo_Symbol** for each gene, each one divided from the next by a space (or a tab).
* Finally there are is a line for each sample that begins with the **Sample_Barcode** and it is followed by the **values** of the given row of the boolean matrix (again spaced by a space or a tab).

## Preprocessor

It is possible to use the **MAF format preprocessor** to convert a standard **MAF** (Multiple annotation Format) file to the **data matrix**; this utility also lets the user to filter the kinds of muytation to be used for the analysis. There are three filtering modes:

* **Standard mode**: this setting automatically selects all the mutations compatible with BML, as described in the publication.
* **All mode**: this setting use all the mutations from the MAF file without any selection. The use of this mode is discouraged.
* **Manual mode**: if none of the check box for the previous modes is checked, the mutations can be manually selected. Write each mutation type using the MAF standard, they can be spaced by any wihespace character.

Similary, it is possible to do the same with the **genes** considered in the analysis.

## BML Execution 

To run **BML** is necessary to compile the form of the main window. We will now describe each field:

* **Job ID**: a mnemonic name for the operation, it will be used in the output files' name, must not be empty.
* **Data Matrix' Path**: the path to the input data matrix, it is possible to look for it by using the 'browse' button. 
* **Number of threads to be used**: it is possible to enable multithreading, this number should be equal to the number of logical CPUs of the machine in use (it is automatically set to this value). It works by parallelizing the tree reseeds and the bootsrap replicates creation and analysis.
* **Cutoff**: minimum number of mutations in samples to include a gene in the analysis. This value can be computed automatically using the convention used in BML publication so that the considered genes are less then the number of samples.
* **Number of tree reseeds**: BML is based on local optimum search. This value indicates the number of times the Bayesan network should be learnt form randomly generated trees. (default is 100)
* **Thresold for inferring paths**: this parameter is used to limit the shown relation in the output tree to the ones with a conditional probability grater than this value. (default is 0.3)
* **Number of bootstrap replicates**: number of times to repeat the bootstrap process mainly to compute the confidence intervals of the probabilities of the relations shown in the output tree. This step is optional and can be disabled. 

> **BML** will outomatically create an **'output'** folder in its directory, the output files will be placed there and have the choosen job ID in their name. 

## Additional informations
To view the output tree it is possible to use the open source **'graphviz'** software (.dot format), however, if an internet connection is available, it can be also viewed by selecting **'graphviz online'** entry in the **'help'** menu and copy and paste the content of the .dot file.