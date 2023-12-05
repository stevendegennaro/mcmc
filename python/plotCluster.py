import mcmc_plot_functions as mc
import structures as st
from parameters import Params
import numpy as np
import pandas as pd
import sys
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RangeSlider,TextBox

############################
##### Create Functions #####
############################

def createIso():
    # Create isochrone with cluster parameters
    st.setClusterParams(pIso,params,pop)
    st.populateIso(pIso)
    global iso
    iso = st.clusterToDataFrame(pIso)
    iso["pop"] = pop
    iso = iso[iso["V"]<90]
    iso = iso[iso["mass1"] >= pIso.contents.evoModels.minMass]
    # with pd.option_context('display.max_rows', None,):
    #   print(iso)

# Create cluster with stars drawn from IMF
def createSim():
    st.setClusterParams(pSim,params,pop)
    pSim.contents.nStars = maxStars
    st.c.reallocStarsArray(pSim)
    st.populateCluster(pSim,params.percentBinary[pop]/100.)
    global sim
    sim = st.clusterToDataFrame(pSim)
    sim["pop"] = pop

def moveSim():
    st.changeClusterParams(pSim,params,pop)
    st.c.evolve(pSim,-1)
    global sim
    sim = st.clusterToDataFrame(pSim)
    sim["pop"] = pop

def createFieldStars():
    minV = min(sim.loc[sim["V"] < 97,"V"])
    maxV = max(sim.loc[sim["V"] < 97,"V"])
    global field
    field = st.populateFieldStars(params,minV,maxV)
    if field is not None:
        field["pop"] = 9

def createScatter():
    if 'sim' not in globals():
        createSim()
    if 'field' not in globals():
        createFieldStars()
    # Scatter the data
    global scattered
    scattered = st.scatterPhot(sim[0:params.nSystems[0]],params)
    if field is not None:
        scatterfield = st.scatterPhot(field,params)
        scattered = pd.concat([scattered,scatterfield])
    #scattered = st.scatterClean(scattered,params)
    #with pd.option_context('display.max_rows', None,): print(scattered.loc[:,"U":"sigI"])

def scatterField():
    global scattered
    # Get rid of the old field stars
    scattered = scattered.loc[scattered["pop"] != 9]

    # Scatter the new field star data only
    if field is not None: 
        scatterfield = st.scatterPhot(field,params)
        scattered = pd.concat([scattered,scatterfield])
    #scattered = st.scatterClean(scattered,params)
    #with pd.option_context('display.max_rows', None,): print(scattered.loc[:,"U":"sigI"])

##########################
#### Change Functions ####
##########################

def newIso(val):
    deleteIso()
    deleteSim()
    deleteScatter()
    deleteFieldStars()
    createIso()
    createSim()
    createFieldStars()
    createScatter()
    plotIso()
    plotSim()
    plotFieldStars()
    plotScatter()

def newSim(val):
    deleteSim()
    deleteScatter()
    createSim()
    createScatter()
    plotSim()
    plotScatter()

def newScatter(val):
    deleteScatter()
    createScatter()
    plotScatter()

def newFieldStars(val):
    deleteFieldStars()
    deleteScatter()
    createFieldStars()
    scatterField()
    plotFieldStars()
    plotScatter()

############################
##### Toggle Functions #####
############################

def hideColor(button):
    button.ax.set_facecolor(hoverColorOn)
    button.color = buttonColorOn
    button.hovercolor = hoverColorOn
    button.label.set_text("Show")

def resetColor(button):
    button.ax.set_facecolor(hoverColor)
    button.color = buttonColor
    button.hovercolor = hoverColor
    button.label.set_text("Hide")

def toggleIso(val):
    [p[0].set_visible(not p[0].get_visible()) for p in isoPoints]
    if isoPoints[0][0].get_visible(): reDrawIso()
    else: hideColor(isoHide)
    fig.canvas.draw()

def toggleSim(val):
    [p.set_visible(not p.get_visible()) for p in simPoints]
    if simPoints[0].get_visible(): reDrawSim()
    else: hideColor(simHide)
    fig.canvas.draw()

def toggleField(val):
    [p.set_visible(not p.get_visible()) for p in fieldPoints]
    if fieldPoints[0].get_visible(): reDrawFieldStars()
    else: hideColor(fieldHide)
    fig.canvas.draw()

def toggleScatter(val):
    scatterPoints[0].set_visible(not scatterPoints[0].get_visible())
    if len(scatterPoints) > 1:
        ebars = scatterPoints[1]
        ebars[0].set_visible(not ebars[0].get_visible())
        for p in ebars[1]:
            p.set_visible(not p.get_visible())
        for p in ebars[2]:
            p.set_visible(not p.get_visible())

    if scatterPoints[0].get_visible(): reDrawScatter()
    else: hideColor(scatterHide)
    fig.canvas.draw()


##############################
##### Plotting functions #####
##############################

def plotBrightLimits():
    global brightLines
    brightLines = [plt.axhline(params.brightLimit[0],linestyle = "dashed",color="gray",zorder=-1)]
    brightLines.append(plt.axhline(params.brightLimit[1],linestyle = "dashed",color="gray",zorder=-1))

def refreshLimits():
    # Figure out min and max x and y values for the graph
    (minxiso,maxxiso),(minyiso, maxyiso) = mc.getHRlims(iso,band1,band2)
    (minxscat,maxxscat),(minyscat, maxyscat) = mc.getHRlims(st.scatterClean(scattered,params),band1,band2,ebars = True)
    if field is not None: (minxfield,maxxfield),(minyfield, maxyfield) = mc.getHRlims(field,band1,band2,ebars = True)
    else: (minxfield,maxxfield),(minyfield, maxyfield) = (99,-99),(-99,99)
    minx = min(minxscat,minxiso,minxfield)
    maxx = max(maxxscat,maxxiso,maxxfield)
    miny = max(minyscat,minyiso,minyfield)
    maxy = min(maxyscat,maxyiso,maxyfield)
    ax.set_xlim([minx - 0.1, maxx + 0.1])
    ax.set_ylim([miny + 1, maxy - 1])
    if 'brightLimit_slider' in globals():
        brightLimit_slider.valmin = maxy
        brightLimit_slider.valmax = miny
        brightLimit_slider.ax.set_xlim(brightLimit_slider.valmax,brightLimit_slider.valmin)

### Isocrhone ###

def plotIso():
    global isoPoints
    if "isoHide" in globals():
        resetColor(isoHide)
    isoPoints = mc.HR(iso, band1 = band1, band2 = band2, type = "l", overplot=True)
    refreshLimits()
    fig.canvas.draw()

def deleteIso():
    [p[0].remove() for p in isoPoints]

# redraws the isochrone without recalculating it
def reDrawIso():
    deleteIso()
    plotIso()

### Simulated stars ###

def plotSim():
    global simPoints
    if "simHide" in globals():
        resetColor(simHide)
    simPoints = mc.HR(sim[0:params.nSystems[0]], 
                        band1 = band1, 
                        band2 = band2, 
                        overplot = True,
                        ebars = False, 
                        color="red", 
                        size = 10,
                        marker = "o")
    refreshLimits()
    fig.canvas.draw()

def deleteSim():
    [p.remove() for p in simPoints]

# Redraws simluated stars without resimulating them
def reDrawSim():
    deleteSim()
    plotSim()

### Scattered stars ###

def plotScatter():
    global scatterPoints
    if "scatterHide" in globals():
        resetColor(scatterHide)
    scatterCleaned = st.scatterClean(scattered,params)
    scatterPoints = mc.HR(scatterCleaned, 
                            band1 = band1, 
                            band2 = band2, 
                            overplot=True, 
                            ebars = True, 
                            size = 10,
                            color = "blue")
    refreshLimits()
    fig.canvas.draw()

def deleteScatter():
    scatterPoints[0].remove()
    if len(scatterPoints) > 1:
        ebars = scatterPoints[1]
        ebars[0].remove()
        for p in ebars[1]:
            p.remove()
        for p in ebars[2]:
            p.remove()

# redraws the scattered points without resimulating or rescattering
def reDrawScatter():
    deleteScatter()
    plotScatter()

### Field Stars ###

def plotFieldStars():
    global fieldPoints
    if "fieldHide" in globals():
        resetColor(fieldHide)
    if field is not None: 
        fieldPoints = mc.HR(field, 
                            band1 = band1, 
                            band2 = band2, 
                            overplot=True, 
                            ebars = True, 
                            color = "green",
                            size= 10,
                            marker = "o")
    refreshLimits()
    fig.canvas.draw()

def deleteFieldStars():
    if field is not None:
        [p.remove() for p in fieldPoints]

# Redraws field stars without resimulating them
def reDrawFieldStars():
    deleteFieldStars()
    plotFieldStars()


#################################
####### Button Functions ########
#################################

def toggleBand1(val):

    global b1,b2,band1

    ### Check whatever band we just clicked. If it's already the band that 
    ### we are using, or the same band as the other band, do nothing
    if band1ax.index(val.inaxes) == b2 or band1ax.index(val.inaxes) == b1:
        return

    ### Otherwise, set band 1 to the band we just clicked
    b1 = band1ax.index(val.inaxes)
    band1 = filters[b1]

    # Set button colors
    for i in range(len(band1ax)):
        if i == b1:
            band1ax[i].set_facecolor(hoverColorOn)
            band1Buttons[i].color = buttonColorOn
            band1Buttons[i].hovercolor = hoverColorOn
        else:
            band1ax[i].set_facecolor(buttonColor)
            band1Buttons[i].color = buttonColor
            band1Buttons[i].hovercolor = hoverColor

    # Update axes and bright/faint limits
    setBands()
    
    # Redraw graphs using new bands
    if isoPoints[0][0].get_visible(): reDrawIso()
    if simPoints[0].get_visible(): reDrawSim()
    if fieldPoints[0].get_visible(): reDrawFieldStars()
    if scatterPoints[0].get_visible(): reDrawScatter()

def toggleBand2(val):

    global b1,b2,band2

    ### Check whatever band we just clicked. If it's already the band that 
    ### we are using, or the same band as the other band, do nothing
    if band2ax.index(val.inaxes) == b2 or band2ax.index(val.inaxes) == b1:
        return

    ### Otherwise, set band 1 to the band we just clicked
    b2 = band2ax.index(val.inaxes)
    band2 = filters[b2]

    # Set button colors
    for i in range(len(band2ax)):
        if i == b2:
            band2ax[i].set_facecolor(hoverColorOn)
            band2Buttons[i].color = buttonColorOn
            band2Buttons[i].hovercolor = hoverColorOn
        else:
            band2ax[i].set_facecolor(buttonColor)
            band2Buttons[i].color = buttonColor
            band2Buttons[i].hovercolor = hoverColor

    # Update axes and bright/faint limits
    setBands()
    
    # Redraw graphs using new bands
    if isoPoints[0][0].get_visible(): reDrawIso()
    if simPoints[0].get_visible(): reDrawSim()
    if fieldPoints[0].get_visible(): reDrawFieldStars()
    if scatterPoints[0].get_visible(): reDrawScatter()

### Function to set band names to the corresponding filter names
### Changes the axes labels and bright/faint limits accordingly
def setBands():
    global b1, b2, band1, band2
    band1 = filters[b1]
    band2 = filters[b2]
    ax.set_ylabel(band2)
    ax.set_xlabel(band1 + " - " + band2)
    params.brightLimit[2] = b2
    fig.canvas.draw()

### Function to set the age and metallicity limits of the sliders
### when the model set is changed
def setModelLimits():
    age_slider.valmin = st.c.getAgeLimit(0)
    age_slider.valmax = st.c.getAgeLimit(1)
    age_slider.ax.set_xlim(age_slider.valmin,age_slider.valmax)
    feh_slider.valmin = st.c.getFeHLimit(0)
    feh_slider.valmax = st.c.getFeHLimit(1)
    feh_slider.ax.set_xlim(feh_slider.valmin,feh_slider.valmax)
    fig.canvas.draw()

### Switch to a different model set
def change_model(val):
    global modelAxes

    ### If we are already using this model set, do nothing
    if modelAxes.index(val.inaxes) == params.msRgbModels:
        return

    # Otherwise, set model set to new model
    modelAxes[params.msRgbModels].set_facecolor(buttonColor)
    modelButtons[params.msRgbModels].color = buttonColor
    modelButtons[params.msRgbModels].hovercolor = hoverColor

    params.msRgbModels = modelAxes.index(val.inaxes)
    st.switchModels(pIso,params)

    #make sure the models are loaded
    fault = st.c.setGlobals(pIso, params.seed, params.verbose)
    assert not fault,f"Problem loading model set {params.msRgbModels}"

    # check to make sure Fe/H and logAge are in range for new models
    if(params.logClusAge[0] < st.c.getAgeLimit(0)):
        params.logClusAge[0] = st.c.getAgeLimit(0)
        age_slider.set_val(params.logClusAge[0])
    if(params.logClusAge[0] > st.c.getAgeLimit(1)):
        params.logClusAge[0] = st.c.getAgeLimit(1)
        age_slider.set_val(params.logClusAge[0])
    if(params.Fe_H[0] < st.c.getFeHLimit(0)):
        params.Fe_H[0] = st.c.getFeHLimit(0)
        feh_slider.set_val(params.Fe_H[0])
    if(params.Fe_H[0] > st.c.getFeHLimit(1)):
        params.Fe_H[0] = st.c.getFeHLimit(1)
        feh_slider.set_val(params.Fe_H[0])
    setModelLimits()

    modelAxes[params.msRgbModels].set_facecolor(hoverColorOn)
    modelButtons[params.msRgbModels].color = buttonColorOn
    modelButtons[params.msRgbModels].hovercolor = hoverColorOn

    st.switchModels(pSim,params)

    deleteIso()
    deleteSim()
    deleteScatter()
    createIso()
    createSim()
    createScatter()
    plotIso()
    plotSim()
    plotScatter()

##########################
#### Slider Functions ####
##########################

def change_nStars(val):
    params.nSystems[0] = int(nStars_slider.val)
    reDrawSim()
    deleteScatter()
    createScatter()
    plotScatter()
    params.updateFileStem()
    save_text_box.set_val(params.file_stem)

def change_nFieldStars(val):
    params.nFieldStars = int(nFieldStars_slider.val)
    newFieldStars(0)
    params.updateFileStem()
    save_text_box.set_val(params.file_stem)

def change_percentBinary(val):
    params.percentBinary[0] = binary_slider.val
    newSim(0)

def change_age(val):
    params.logClusAge[0] = age_slider.val
    newIso(0)
    moveSim()
    reDrawSim()
    refreshLimits()
    params.updateFileStem()
    save_text_box.set_val(params.file_stem)

def change_feh(val):
    params.Fe_H[0] = feh_slider.val
    newIso(0)
    moveSim()
    reDrawSim()
    refreshLimits()
    params.updateFileStem()
    save_text_box.set_val(params.file_stem)

def change_mod(val):
    params.distMod[0] = mod_slider.val
    deleteIso()
    createIso()
    plotIso()
    moveSim()
    reDrawSim()
    deleteScatter()
    createScatter()
    plotScatter()

def change_av(val):
    params.Av[0] = av_slider.val
    deleteIso()
    createIso()
    plotIso()
    moveSim()
    reDrawSim()
    deleteScatter()
    createScatter()
    plotScatter()

def change_exposure(new_val,i):
    params.exposures[i] = new_val
    deleteScatter()
    createScatter()
    plotScatter()

def change_s2n(new_val):
    params.limitS2N = new_val
    reDrawScatter()

def change_brightLimit(new_val):
    for line in brightLines:
        line.remove()
    params.brightLimit[0] = new_val[0]
    params.brightLimit[1] = new_val[1]
    plotBrightLimits()
    reDrawScatter()

##############
#### Main ####
##############

pop = 0
maxStars = 1000

# Initialize
params = Params()
np.random.seed(params.seed)
pIso = st.createCluster(params)
pSim = st.createCluster(params)
pField = st.createCluster(params)
fault = st.c.setGlobals(pIso, params.seed, params.verbose)
if fault: sys.exit()
st.setClusterParams(pIso,params,pop)
createIso()
createSim()
createFieldStars()
createScatter()

# Set labels and make room for buttons and sliders
plt.ion()
fig,ax = plt.subplots(figsize=(13,7))
figwidth = 0.45
fig.subplots_adjust(right=figwidth,left=0.05,top=0.95)
filters = st.getFilters()
b1 = 1; b2 = 2
setBands()
plotIso()
plotSim()
plotScatter()
plotBrightLimits()
plotFieldStars()

# Button sizes
buttonheight = 0.05
buttonwidth = 0.1
pad = 0.01

# Add filter buttons
filterWidth = (2 * buttonwidth + pad) / 8.
band1Buttons = [None] * len(filters)
band2Buttons = [None] * len(filters)
band1ax      = [None] * len(filters)
band2ax      = [None] * len(filters)
for i,filt in enumerate(filters):
    band1ax[i] = plt.axes([figwidth + pad + i * filterWidth, .95 - (buttonheight), filterWidth, buttonheight])
    band1Buttons[i] = Button(band1ax[i], filt)
    band1Buttons[i].on_clicked(toggleBand1)
    band2ax[i] = plt.axes([figwidth + pad + i * filterWidth, .95 - (2 * buttonheight), filterWidth, buttonheight])
    band2Buttons[i] = Button(band2ax[i], filt)
    band2Buttons[i].on_clicked(toggleBand2)

# style filter buttons
buttonColor = band1Buttons[0].color
hoverColor = band1Buttons[0].hovercolor
buttonColorOn = "dimgray"
hoverColorOn = "darkgray"
band1ax[b1].set_facecolor(buttonColorOn)
band1Buttons[b1].color = buttonColorOn
band1Buttons[b1].hovercolor = hoverColorOn
band2ax[b2].set_facecolor(buttonColorOn)
band2Buttons[b2].color = buttonColorOn
band2Buttons[b2].hovercolor = hoverColorOn

#### Buttons to Regenerate and to Hide/show
# Isochrones
axes = plt.axes([figwidth + pad, .95 - (3 * buttonheight + pad), buttonwidth, buttonheight])
isoButton = Button(axes, 'Isochrone')
isoButton.on_clicked(newIso)
axes = plt.axes([figwidth + pad + buttonwidth + pad, .95 - (3 * buttonheight + pad), buttonwidth, buttonheight])
isoHide = Button(axes, 'Hide')
isoHide.on_clicked(toggleIso)

# sim
axes = plt.axes([figwidth + pad, .95 - (4 * buttonheight + 2 * pad), buttonwidth, buttonheight])
simButton = Button(axes, 'simulate')
simButton.on_clicked(newSim)
axes = plt.axes([figwidth + pad + buttonwidth + pad, .95 - (4 * buttonheight + 2 * pad), buttonwidth, buttonheight])
simHide = Button(axes, 'Hide')
simHide.on_clicked(toggleSim)

# field stars
axes = plt.axes([figwidth + pad, .95 - (5 * buttonheight + 3 * pad), buttonwidth, buttonheight])
fieldButton = Button(axes, 'field stars')
fieldButton.on_clicked(newFieldStars)
axes = plt.axes([figwidth + pad + buttonwidth + pad, .95 - (5 * buttonheight + 3 * pad), buttonwidth, buttonheight])
fieldHide = Button(axes, 'Hide')
fieldHide.on_clicked(toggleField)

# scatter
axes = plt.axes([figwidth + pad, .95 - (6 * buttonheight + 4 * pad), buttonwidth, buttonheight])
scatterButton = Button(axes, 'scatter')
scatterButton.on_clicked(newScatter)
axes = plt.axes([figwidth + pad + buttonwidth + pad, .95 - (6 * buttonheight + 4 * pad), buttonwidth, buttonheight])
scatterHide = Button(axes, 'Hide')
scatterHide.on_clicked(toggleScatter)

#### Model Switching Buttons ####
modelNames = ['Girardi','Chaboyer','Yale-Yonsai','DSED']
modelPosition = [[0,0],[1,1],[0,1],[1,0]]
modelButtons = [None] * 4
modelAxes = [None] * 4

for i in range(4):
    modelAxes[i] = plt.axes([figwidth + pad + modelPosition[i][0] * (buttonwidth + pad), 
                             .95 - (10 * buttonheight + pad + modelPosition[i][1] * (buttonheight + pad)), 
                             buttonwidth, 
                             buttonheight])
    modelButtons[i] = Button(modelAxes[i], modelNames[i])
    modelButtons[i].on_clicked(change_model)

modelAxes[params.msRgbModels].set_facecolor(buttonColorOn)
modelButtons[params.msRgbModels].color = buttonColorOn
modelButtons[params.msRgbModels].hovercolor = hoverColorOn


#################
#### Sliders ####
#################
slider_width = (1.8 * buttonwidth)
slider_height = 0.04

# N Stars
param_ax = fig.add_axes([figwidth + (2 * buttonwidth) + (8 * pad), .95 - (slider_height), slider_width, slider_height])
nStars_slider = Slider(
    ax=param_ax,
    label="N",
    valmin=1,
    valmax=maxStars,
    valinit=params.nSystems[0],
    valstep = 1
)
nStars_slider.on_changed(change_nStars)

# N field Stars
param_ax = fig.add_axes([figwidth + (2 * buttonwidth) + (8 * pad), .95 - (2 * slider_height), slider_width, slider_height])
nFieldStars_slider = Slider(
    ax=param_ax,
    label="FS",
    valmin=0,
    valmax=1000,
    valinit=params.nFieldStars,
)
nFieldStars_slider.on_changed(change_nFieldStars)

# percent Binary
param_ax = fig.add_axes([figwidth + (2 * buttonwidth) + (8 * pad), .95 - (3 * slider_height), slider_width, slider_height])
binary_slider = Slider(
    ax=param_ax,
    label="% Binary",
    valmin=0,
    valmax=100,
    valinit=params.percentBinary[0],
)
binary_slider.on_changed(change_percentBinary)

# Log Age
param_ax = fig.add_axes([figwidth + (2 * buttonwidth) + (8 * pad), .95 - (4 * slider_height), slider_width, slider_height])
age_slider = Slider(
    ax=param_ax,
    label="Log Age",
    valmin=0,
    valmax=10,
    valinit=params.logClusAge[0],
)
age_slider.on_changed(change_age)

# FeH
param_ax = fig.add_axes([figwidth + (2 * buttonwidth) + (8 * pad), .95 - (5 * slider_height), slider_width, slider_height])
feh_slider = Slider(
    ax=param_ax,
    label="[Fe/H]",
    valmin=0,
    valmax=10,
    valinit=params.Fe_H[0],
)
feh_slider.on_changed(change_feh)

setModelLimits()

# modulus
param_ax = fig.add_axes([figwidth + (2 * buttonwidth) + (8 * pad), .95 - (6 * slider_height), slider_width, slider_height])
mod_slider = Slider(
    ax=param_ax,
    label="m-M",
    valmin=0,
    valmax=20,
    valinit=params.distMod[0],
)
mod_slider.on_changed(change_mod)

# Absorption
param_ax = fig.add_axes([figwidth + (2 * buttonwidth) + (8 * pad), .95 - (7 * slider_height), slider_width, slider_height])
av_slider = Slider(
    ax=param_ax,
    label="Av",
    valmin=0,
    valmax=1.0,
    valinit=params.Av[0],
)
av_slider.on_changed(change_av)

# Exposures
exposure_sliders = []
exposure_axes = [None] * len(filters)
for i,filt in enumerate(filters):
    # Absorption
    exposure_ax = fig.add_axes([figwidth + (2 * buttonwidth) + (8 * pad), .95 - ((9 + i) * slider_height), slider_width, slider_height])
    exposure_sliders.append(
        Slider(
            ax=exposure_ax,
            label=filt,
            valmin=0,
            valmax=0.2,
            valinit=params.exposures[i],
        )
    )
for s in range(len(exposure_sliders)):
    exec(f"exposure_sliders[{s}].on_changed(lambda new_val: change_exposure(new_val, {s}))")

# Limiting Signal to Noise
s2n_ax = fig.add_axes([figwidth + (2 * buttonwidth) + (8 * pad), .95 - ((17) * slider_height), slider_width, slider_height])
s2n_slider = Slider(
    ax=s2n_ax,
    label="limitS2N",
    valmin=0,
    valmax=20,
    valinit=params.limitS2N,
    closedmin = False
)
s2n_slider.on_changed(change_s2n)

# Bright and Faint Limits
brightLimit_ax = fig.add_axes([figwidth + (2 * buttonwidth) + (8 * pad), .95 - ((18) * slider_height), slider_width, slider_height])
brightLimit_slider = RangeSlider(
    ax=brightLimit_ax,
    label="limits",
    valmin=0,
    valmax=30,
    valinit=(params.brightLimit[0],params.brightLimit[1]),
    valfmt = "%.1f"
)
brightLimit_slider.on_changed(change_brightLimit)
refreshLimits()

#####################
#### File Saving ####
#####################

def getScatterHeader():
    header = ["id"]
    for filt in range(st.FILTS):
        if st.c.getUseFilt(filt):
            if params.exposures[filt] > 0.0:
                header.append(f"{st.c.getFilterName(filt).decode()}")
    for filt in range(st.FILTS):
        if st.c.getUseFilt(filt):
            if params.exposures[filt] > 0.0:
                header.append(f"sig{st.c.getFilterName(filt).decode()}")
    header.extend(["mass1", "massRatio", "stage1", "CMprior", "useDuringBurnIn", "pop"])
    return header

def save_sim(val):
    filename = save_text_box.text + ".out"
    sim_out = pd.concat([sim[0:params.nSystems[0]],field])

    with open(filename, "w") as f:
        sim_out.to_csv(f,sep=" ",index=False)

def save_scatter(val):
    filename = save_text_box.text + ".scatter"
    header = getScatterHeader()
    scatter_out  = st.scatterClean(scattered,params)

    scatter_out["massRatio"] = scatter_out["mass2"]/scatter_out["mass1"]
    clusterMemberPrior = len(scatter_out[scatter_out["pop"] != 9].index)/len(scatter_out.index)
    scatter_out["CMprior"] = clusterMemberPrior
    scatter_out["useDuringBurnIn"] = 1
    scatter_out.loc[scatter_out["pop"] == 9,"useDuringBurnIn"] = 0
    scatter_out = scatter_out.reindex(columns=header)
    with open(filename, "w") as f:
        scatter_out.to_csv(f,sep=" ",index=False)

textwidth = 1 - (figwidth + 3 * pad + 2 * buttonwidth) - pad
axes = plt.axes([1 - textwidth - pad, ax.get_position().bounds[1], textwidth, buttonheight])
save_text_box = TextBox(axes,label = None,initial = params.file_stem)

axes = plt.axes([figwidth + pad, ax.get_position().bounds[1], buttonwidth, buttonheight])
save_sim_button = Button(axes, 'Save Sim')
axes = plt.axes([figwidth + pad + buttonwidth + pad, ax.get_position().bounds[1], buttonwidth, buttonheight])
save_scatter_button = Button(axes, 'Save Scatter')
save_sim_button.on_clicked(save_sim)
save_scatter_button.on_clicked(save_scatter)

plt.sca(ax)
plt.show()



