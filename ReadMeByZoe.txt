___________________________________ Instructions to use the codes for the EEC analysis ___________________________________

First of all, we need to understand what we want to analyze. 
We can want the EEC(dr) distribution or E3C(dr) one. We will focus on the EEC for now.

To get the EEC distribution we need to follow these main 3 steps:
1) Get the EEC histograms
2) Perform the Template Fit
3) Unfold data

and voilà.

We can also do a lot of plots to understand the influence of different variables in the EEC(dr). We'll see them when needed.

_____________________________________________________________________________________________________________________________

1) GET EEC 3D (get_eec_3d.cpp)

What does this code? It generates the eec distributions at reco level (for data & MC) and gen level (only for MC). Also, for MC, the efficiency histograms and tag efficiency histograms.

This code needs three things to run:
a) The "tTree.h"
b) The "binning_hist.h"
c) Lida's data (path in the code)

Regarding a), I modified it from Francesca to use the same .h for MC/data. Also, I changed b) to solve some bugs. You can change the number of bins there manually. And maybe it would be a good idea to do a backup of Lida's data in someone's directory.

Now, we can select what we want to run choosing the following parameters:

- dataType: -1 for LowEG data, 0 for HighEG data, 1 for bjet, 2 for dijet. It will be display in the name of the histos as LowEG/HighEG_data or bjet/dijet_MC.
- btag: true/false. In what Francesca did was almost always true.
- aggregation: true/false.
- ideal aggregation: if aggregation is true, it can be true/false.
- pT range: in some parts you set the lower and upper limit, by default 80_140.
- folder: we put MANUALLY the folder where we want to save the histograms.

The code will run for all flavours in case it's MC, even though this is not applicable for bjets. 
The histograms will be empty if it doesn't apply.

So, you will have one file for each combination of parameters with all the histograms inside.

For data/all the MC, the histogram is called "h_dataType_all". For the flavours, "h_dataType_flavour", for instance "h_MC_dijet_light". Our flavours are "all", "b1" (only 1 B hadron), "b2" (two or more B hadrons, called "moreb" by Francesca), "c" and "light" (which's sum was called by Francesca "other").

In the last lines of the code, you can comment if you don't one one or more of the histograms. The efficiencies one are not optimized.
To check the progress, it should print the number of events already analyzed 10 times, so each one is 10%.

Moreover, you can also do the E3C histograms with this code and that part is also optimized (you have one function for the reco and other for the gen).

EXTRA: (or not so extra)
If you change only one parameter at a time, you can see its influence. The variations Francesca studied are
- using n = 1 or n = 2.
- effect on pT. 
- btag vs no btag (should be made with ideal aggregation to have the same as Francesca I think.)
- flavour effect (for dijet, no b tagged, aggregated)
and more...

These plots can also be used as control plots, in the sense that we can check if different aspects of the data are the same as Francesca's.

_____________________________________________________________________________________________________________________________

2) DO TEMPLATE FIT 3D (do_template_fit_3d.cpp)

In this case, you need to include:
a) "tTree.h"
b) "TFile.h" (if not, it crashes)

The main functions of the code are:
- do_template_fit: calculates template fit
- draw_template_fit_result: draws it
- draw_eec: draws the EEC distribution
- do_template_fit_3d: the macro, compiles previous functions. We put as input if it's MC or not (default: isMC = false).

This code is much simpler but you have to take into account some important things:
- The binning you chose before will be VERY important. For small bins (less than dr = 0.1), the fit would not converge. So that's why I have data regenerated for "smaller bins". I have less binning region and also less number of bins in the region I'm considering.
- The macro function can be optimized more.
- All the functions have been optimized for the new data formatting. If you try to run it with Francesca's data, it won't work since she uses one file per histogram and here we use one file for each dataType.

_____________________________________________________________________________________________________________________________

3) DO THE UNFOLDING (you will need several codes for this, so it will be divided in smaller parts)


