# EAO code
This is The repository for the code in the paper  about approximating overlapping group lasso by non-overlapping group lasso.

There are two folders in this file, "simulation" and "application", corresponding to the simulation examples and real data analysis in the paper.

The folder "simulation" contains two sub-folders, "interlocking" and "gene pathway". Each sub-folder is for a specific group strucuture setting in the paper. Runing through the main functions in each sub-folder will generate the results in the simulation section. All of the 
auxiliary functions for the main function are in the "sub-function" folder. Users need to add approprate path in matlab to merge these functions together.

The folder "application" is the code for real data analysis. Similary, runing the main functions generate the results in the paper, and all of the auxillary functions are in the "sub-function" folder.

Additionally, our code also used some functions in the "SLEP" package (see J. Liu, S. Ji, and J. Ye. SLEP: Sparse Learning with Efficient Projections.).
The "SLEP"  package could be downloaded here: http://yelabs.net/software/SLEP/.  Users also need to add approprate path in matlab to merge these functions together.

The code was originally generately by the first author of the paper, Mingyu Qi. Please contact us for questions (mq3sq@virginia.edu, tianxili@virginia.edu).
