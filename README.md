# Fifty shades of grey

This repository presents the implementation section of the journal article [Fifty shades of grey: Automated stochastic model identification of building heat dynamics](https://doi.org/10.1016/j.enbuild.2021.111195) examplified on a toy-set of 3 anonymized buildings.


## Table of content 

-   [Directory structure](#directory-structure)
-   [Model identification](#model-identification)
    -   [Model selection](#model-selection)
    -   [Model evaluation](#model-evaluation)
-   [Author](#authors)
-   [Licence](#license)


## Directory structure

The directory is organized given the following structure;
-   [src](https://github.com/JulienLeprince/fiftyshadesofgrey/tree/main/src)
    -   [RC models](https://github.com/JulienLeprince/fiftyshadesofgrey/blob/main/src/allmodels.R)
    -   [Model selection](https://github.com/JulienLeprince/fiftyshadesofgrey/blob/main/src/main1_modelselection_parallel.R)
    -   [Results processing](https://github.com/JulienLeprince/fiftyshadesofgrey/blob/main/src/main2_resultsprocessing.R)
    -   [Model evaluation](https://github.com/JulienLeprince/fiftyshadesofgrey/blob/main/src/main3_modelevaluation.R)
    -   [Utility functions](https://github.com/JulienLeprince/fiftyshadesofgrey/blob/main/src/utils.R)
-   [data](https://github.com/JulienLeprince/fiftyshadesofgrey/tree/main/src/data)
    -   [in](https://github.com/JulienLeprince/fiftyshadesofgrey/tree/main/src/data/in)
    -   [out](https://github.com/JulienLeprince/fiftyshadesofgrey/tree/main/src/data/out)
-   [figures](https://github.com/JulienLeprince/fiftyshadesofgrey/tree/main/fig)


To get started, simply clone this repository on your computer or Fork it via GitHub. After installing dependencies from  the `requirements.txt` file, the code should run properly.


## Model selection

This study considers grey-box models ranging from a simple first order model \textit{Ti}, where the inside temperature state-point $T_i$ and its RC parameters $R_{ia}$ and $C_i$ are solely treated, to 5$^{th}$ order ones, where the addition of sensor $T_s$, medium $T_m$, heater $T_h$ and building envelope $T_e$ state points along with their respective RC parameters each add a variety of model extensions to chose from. Additionally, the building envelope component proposes additional parameter extensions modeling direct inside to outside heat exchanges and facade solar gains, which are here considered as a block extension $A_eR_{ia}$.
For a detailed description of the grey-box models, the reader is suggested to refer to the work of [Bacher and Madsen](https://www.sciencedirect.com/science/article/pii/S0378778811000491).

The model selection procedure employs a likelihood ratio test to statistically determine whether a more complex model performs significantly better, or not, compared to a simpler, sub-model.
A forward selection procedure is proposed beginning with the simplest feasible model, $T_i$, and extending it iteratively with the component presenting the most significant improvement. The procedure terminates when no model extension yields a p-value below the pre-specified limit, commonly fixed at 5\%.
Possible candidates for model improvement are selected from a set of predefined extensions, resulting from the combination possibilities of the different considered model components, i.e., $T_e$, $T_h$, $T_m$, $T_s$, $A_eR_{ia}$. 
The overall model selection scheme is presented below. Possible model combinations are mapped and linked, visually exposing the different existing paths of the forward selection procedure.

<p align="center">
  <img src="https://github.com/JulienLeprince/fiftyshadesofgrey/tree/main/fig/modelselection.png" alt="Model selection"/>
</p>

Grey-box models are here implemented using the computer software [CTSM-R](http://ctsm.info/) developed at the Technical University of Denmark. It produces maximum likelihood estimates of model parameters thanks to an optimization algorithm performed over a Kalman filter.


### Model evaluation

Models are finally evaluated leveraging the commonly employed qualitative appreciation of model fits from cumulated periodograms of the residuals.
Typically, an appropriate model will produce residuals with Gaussian white-noise properties, which in the frequency domain, denotes a theoretical constant periodogram. Observing whether obtained model residuals are located around this straight line, e.g. within a surrounding confidence interval, consequently serves as an appropriate indicator of a model's quality.

By calculating the frequency differences between a selected model's Cumulated Periodogram (CP) and its confidence interval, we obtain boundary excess values which, in turn, can be summed into a unique numerical indicator, i.e., the Cumulated Periodogram Boundary Excess Sum (CPBES). This indicator characterizes the amount of auto-correlation present in the considered time-series, which implies white noise properties when close to zero. CPBES consequently allows the differentiation of poor, suitable and good models resulting from the previous forward selection procedure. To allow fair comparisons of CPBES between times series of different length, we normalize it by its length, and obtain the normalized CPBES (nCPBES).

<p align="center">
  <img src="https://github.com/JulienLeprince/fiftyshadesofgrey/tree/main/fig/nCPBES_demo_final.png" alt="nCPBES Demo"/>
</p>


## Authors

**[Julien Leprince](https://github.com/JulienLeprince)**

**Prof. [Henrik Madsen](https://henrikmadsen.org/)**

**[Jaume Palmer Real](https://orbit.dtu.dk/en/persons/jaume-palmer-real)**

**Prof. [Clayton Miller](https://github.com/cmiller8)**

**[Rik van der Vlist](https://www.linkedin.com/in/rik-van-der-vlist-124b62138/)**

**Dr. [Kaustav Basu](https://www.linkedin.com/in/kaustav-basu-phd-5973311b/)**

**Prof. [Wim Zeiler](https://www.tue.nl/en/research/researchers/wim-zeiler/)**


## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details