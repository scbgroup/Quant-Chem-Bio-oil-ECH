from pylab import *
from energydiagram import ED

'''
These are the 2 functions needed to produce energy diagrams
'build_diagram' builds the diagram using the ED package
'fig_settings' applies customized features and matplotlib formatting to the figure and attaches the diagram
'''

def build_diagram(rxnE, states, colors, ax):
    dia = ED()
    '''
    # plot the energy of each state on the surfaces
    if len(shape(rxnE)) == 1:
        rxn_l = len(rxnE)
        for i, a in enumerate(array(rxnE).T):
            dia.add_level(a, '', top_text = '', color = colors)
        # add links between states on a specific surface
        
        t0 = -1
        t1 = len(rxnE) - 1
        max_val = rxn_l - 1
        print(max_val)
        for i,ii in enumerate(rxnE):
            for j,jj in enumerate(states):
                if t1 == max_val:
                    break
                t0 += 1
                t1 += 1 
                print(t0,t1)
                dia.add_link(t0,t1)
    
    else:
    '''
    for i, a in enumerate(array(rxnE).T):
        for j, e in enumerate(a):
            # plot the same state in line with other surfaces
            #dia.dimension = dim
            if j > 0:
                dia.add_level(e, '', 'last', top_text = '', color = colors[j])
            else:
                dia.add_level(e, '', top_text = '', color = colors[j])
        
    
    # add links between states on a specific surface
    t0 = -1
    t1 = len(rxnE) - 1
    max_val = (shape(rxnE)[0]*shape(rxnE)[1])-1
    for i,ii in enumerate(rxnE):
        for j,jj in enumerate(states):
            if t1 == max_val:
                break
            t0 += 1
            t1 += 1 
            dia.add_link(t0,t1)#, color = colors[j])
    # the diagram has to be plotted before matplotlib settings can be adjusted
    #dia.plot(ylabel = 'Relative Energy / eV')
    #dia.ratio = 3
    dia.dimension = 3
    dia.space = 1
    #dia.bottom_text_fontsize = 8.5
    #dia.dimension = 50
    #dia.offset = 0.1
    #dia.aspect = 'equal'
    
    dia.plot(ylabel = 'Relative Free Energy / eV', ax = ax)
    return dia

def fig_settings(diagram, ylim, labels, title, st_lab, lb_fnt, x, y):
    # update fig details
    
    # change y axis label
    diagram.ax.set_ylabel('Relative Free Energy / eV')
    
    # set y lim to remove conflicts with labels
    diagram.ax.set_ylim(ylim)
    
    # show x axis and set x axis label
    diagram.ax.axes.get_xaxis().set_visible(True)
    diagram.ax.spines['bottom'].set_visible(True)
    diagram.ax.spines['top'].set_visible(True)
    diagram.ax.spines['right'].set_visible(True)
    diagram.ax.set_xticks([ ])
    diagram.ax.set_xlabel("Reaction Coordinate")
    # diagram.ax.set_title(title)
    # diagram.fig.set_figwidth(8)
    diagram.ax.text(0.5,0.95,title, horizontalalignment='center', verticalalignment='top', transform=diagram.ax.transAxes)

    # set manual state labels
    for i, ii in enumerate(st_lab):
        diagram.ax.text(x[i],y[i], ii,rotation = 90, transform=diagram.ax.transAxes, fontsize = lb_fnt)
    
    
    # add a legend
    # diagram.ax.legend(labels, loc = 'upper right')#bbox_to_anchor=(1.04, 0.5), loc="center left")
    return diagram.fig