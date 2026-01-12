# Developing a Machine Learning Tool for Early Prediction of Lung Cancer Using miRNA Data.

MicroRNAs (miRNAs) are a class of small, non-coding RNA molecules that play a crucial role in the regulation of gene expression at the post-transcriptional level (O'Brien et al., 2018). They exert their regulatory function by binding to the 3' untranslated region (UTR) of target mRNAs modulating a variety of cellular processes, including proliferation, differentiation, and apoptosis, thereby playing a key role in maintaining cellular homeostasis (O'Brien et al., 2018). 
miRNAs are considered a promising biomarker for cancer diagnosis because of their stability in bodily fluids making them a suitable tool for non-invasive diagnosis (Lu et al., 2005; Mitchell et al., 2008). Recently, machine learning (ML) has been used as a powerful tool for analyzing and processing miRNA expression profiles helping in identifying various cancerous cases and cancer stages, including lung cancer. Further, the integration between miRNA data and ML algorithms has played a significant role in overcoming the limitations of traditional diagnostic methods such as low sensitivity and specificity.
The current study explored multiple machine-learning approaches (SVM, RF, GLM) while optimizing feature selection and hyperparameters to identify the most effective miRNA biomarkers for early detection of lung cancer. 

miRNA-based lung cancer prediction workflow included 5 steps 1-Data Acquisition from GEO 2-Data Processing (filtering, normalization and splitting training (80%) and testing (20%)) 3-Feature Selection and Model Development (RF, SVM and GLR) 4-Model Evaluation and Validation (Cross-Validation Framework and ROC and AUCs Metrics) 5-Biological Interpretation (differential expression analysis using Limma).


<img width="448" height="442" alt="image" src="https://github.com/user-attachments/assets/115199fd-7590-41e1-ae1b-ab6d8cc14dcf" />

Figure 1: 3D Principal Component Analysis (PCA) of miRNA expression before and after preprocessing. Left panel shows the original dataset with samples projected along the first three principal components (PC1, PC2, and PC3), colored by group: Lung cancer (red) and Noncancer (blue). Right panel displays the same dataset after low-expression filtering and quantile normalization, demonstrating improved cluster separation between Lung cancer (red) and Noncancer (blue) groups along the three principal components.


<img width="975" height="450" alt="image" src="https://github.com/user-attachments/assets/90ae4333-a856-420d-b04d-7df549982de5" />

Figure 2: Receiver Operating Characteristic (ROC) curves comparing three machine learning models: Support Vector Machine (SVM), Logistic Regression (GLM), and Random Forest (RF) for predicting lung cancer vs. non-cancer based on miRNA expression. Each panel includes the AUC, sensitivity, and specificity at the optimal threshold.


<img width="959" height="581" alt="image" src="https://github.com/user-attachments/assets/0888941c-9787-47ef-bae6-1d2548f61a5e" />

Figure 3: PCA plot using the top miRNAs selected by the Random Forest model, showing the projection of samples onto the first two principal components (PC1 and PC2). Each point represents a sample colored by cancer status, with 95% confidence ellipses overlaid for each group.


<img width="975" height="747" alt="image" src="https://github.com/user-attachments/assets/89245b03-cd9d-48c8-bc62-8d05ab9985d7" />

Figure 4: Calibration plot of the Random Forest model showing the relationship between predicted probabilities and actual observed frequencies for lung cancer classification. The dashed diagonal line represents perfect calibration.


<img width="975" height="745" alt="image" src="https://github.com/user-attachments/assets/4757c29d-0562-437c-b63c-09f809794fb1" />

Figure 5: Barplot showing the top 30 most important miRNAs ranked by their variable importance scores derived from the Random Forest model. Importance scores reflect the contribution of each miRNA to the model’s classification accuracy.


<img width="975" height="597" alt="image" src="https://github.com/user-attachments/assets/629569cb-29b8-48d8-8923-c5eb3a7dbcc9" />

Figure 6: Boxplot showing the normalized expression (log2) of the top 30 most important miRNAs selected by the Random Forest model, comparing expression levels between lung cancer and non-cancer samples.


<img width="975" height="542" alt="image" src="https://github.com/user-attachments/assets/8a332bb9-4283-421d-9c82-75538986b465" />

Figure 7: Volcano plot of differentially expressed miRNAs comparing lung cancer versus non-cancer samples. Top 30 miRNAs selected by Random Forest are annotated. miRNAs are classified as upregulated (red), downregulated (blue), or not significant (gray) based on log2 fold change > 1 and adjusted p-value < 0.05.
