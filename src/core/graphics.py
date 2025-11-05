import matplotlib.pyplot as plt

def blockiness_calculator_graphics(blockiness_calculator):
    """
    Generate visualization of charge blockiness analysis results.
    
    Creates a 2-panel figure showing the raw charge profile and processed
    blockiness calculations for protein sequence analysis.
    
    Args:
        blockiness_calculator: object containing blockiness calculation results
    
    Returns:
        matplotlib.figure.Figure
    """
    bc = blockiness_calculator
    fig, ax = plt.subplots(num = bc.uid, nrows = 2, ncols = 1, figsize=(12,7))

    fig.subplots_adjust(left=.1, right=.9, top=.91, bottom=0.11, hspace=.320)

    ax[0].plot(range(len(bc.seq_charge)), bc.seq_charge, linewidth=.2)

    ax[0].set_xlabel('Residues', fontsize=14)
    ax[0].set_ylabel('Charge ', fontsize=14)
    ax[0].set_title('Charge profile')

    ax[1].plot(range(len(bc.lpf_seq)), bc.lpf_seq, color='darkslateblue', linewidth=.5)
    ax[1].plot(range(len(bc.lpf_lcr)), bc.lpf_lcr, color='plum', alpha=.5)
    ax[1].plot(range(len(bc.blockiness)), bc.blockiness)

    ax[1].set_xlabel('Residues', fontsize=14)
    ax[1].set_ylabel('Charge ', fontsize=14)
    ax[1].set_title('Charge Blockiness')

    return fig