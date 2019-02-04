import  matplotlib;                matplotlib.use('PDF')


def latexify(fig_width=None, fig_height=None, columns=1, equal=False, fontsize=10, ratio=None, ggplot=True, usetex=True):
    '''                                                                                                                                                              
    Set up matplotlib's RC params for LaTeX plotting.                                                                                                                 
    Call this before plotting a figure.                                                                                                                                                                                                                                                                                                    
    Parameters                                                                                                                                                        
    ----------                                                                                                                                                         
    fig_width:  float, optional, inches                                                                                                                               
    fig_height: float, optional, inches                                                                                                                                                                                                                                                                                                     
    Columns : {1, 2}                                                                                                                                                  
    '''

    ##  Code adapted from http://www.scipy.org/Cookbook/Matplotlib/LaTeX_Examples                                                                                     
    ##  Width and max height in inches for IEEE journals taken from                                                                                                   
    ##  computer.org/cms/Computer.org/Journal%20templates/transactions_art_guide.pdf                                                                                   
    import  matplotlib.pyplot as      plt
    import  matplotlib        as      mpl

    from    numpy             import  sqrt

    assert(columns in [1,2])

    if ggplot:
      plt.style.use('ggplot')

    mpl.rc('text', usetex = usetex)

    if fig_width is None:
        fig_width = 6.08948 if columns == 2 else 6.08948 / 2.  ##  width in inches
        
    if fig_height is None:
      if equal is False:
          if ratio is None:
            ratio     = (sqrt(5.) - 1.0) / 2.0

          fig_height  = fig_width * ratio                        ## Height in inches                                                                                     

      else:
        fig_height  = fig_width

    MAX_HEIGHT_INCHES = 12.1

    if fig_height > MAX_HEIGHT_INCHES:
        print("WARNING: fig_height too large:" + str(fig_height) + "so will reduce to" + str(MAX_HEIGHT_INCHES) + "inches.")

        fig_height = MAX_HEIGHT_INCHES

    params   = {'backend': 'pdf',
                ## 'text.latex.preamble': ['\usepackage{gensymb}'],                                                                                                   
                'axes.labelsize': fontsize,  ## fontsize for x and y labels (was 10)                                                                                  
                'axes.titlesize': fontsize,
                'font.size': fontsize,
                'legend.fontsize': fontsize,
                'xtick.labelsize': fontsize,
                'ytick.labelsize': fontsize,
                'text.usetex': usetex,
                'figure.figsize': [fig_width, fig_height],
                'font.family': 'serif',
                'figure.facecolor': 'w'}

    matplotlib.rcParams.update(params)


def sci_notation(num, decimal_digits=1, precision=None, exponent=None):
    '''                                                                                                                                                                   Returns a string representation of the scientific                                                                                                                 
    notation of the given number formatted for use with                                                                                                               
    LaTeX or Mathtext, with specified number of significant                                                                                                           
    decimal digits and precision (number of decimal digits                                                                                                            
    to show). The exponent to be used can also be specified                                                                                                           
    explicitly.                                                                                                                                                                                                                                                                                                                             
    https://stackoverflow.com/questions/18311909/how-do-i-annotate-with-power-of-ten-formatting                                                                       
    '''
    if not exponent:
        exponent = int(floor(log10(abs(num))))

    coeff = round(num / float(10**exponent), decimal_digits)

    if not precision:
        precision = decimal_digits

    return  r"${0:.{2}f}\cdot10^{{{1:d}}}$".format(coeff, exponent, precision)
