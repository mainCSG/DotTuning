import numpy as np
import pandas as pd
import wx

# TODO: input() does not work in Python 2

"""
Get user input for doing virtual sweep
"""

# WX menu to open file dialog for csv
def get_path(wildcard):
    app = wx.App(None)
    style = wx.FD_OPEN | wx.FD_FILE_MUST_EXIST
    dialog = wx.FileDialog(None, 'Open', wildcard=wildcard, style=style)
    if dialog.ShowModal() == wx.ID_OK:
        path = dialog.GetPath()
    else:
        path = None
    dialog.Destroy()
    return path

# Splash
splash = '*******************************************************************************\n'
splash += 'Virtual Gating Setup Script\n'
splash += '---------------------------\n'
splash += 'This program helps setup a virtual sweep to ensure the same region of charge\n'
splash += 'stability diagram is swept.\n'
splash += '*******************************************************************************\n'
print(splash)

# Load results, capacitance matrix
raw_input('Press Enter to load CSV')

# filename = str(raw_input("Capacitance value CSV: "))
filename = get_path('*.csv')

data = pd.read_csv(filename)

print('\nSelect what virtual gating mode to apply:')

skewLabel = '\n[skew]  :  Only gate terms are known: Cg1d1,Cg2d2,Cg1d2,Cg2d1'
fullLabel = '\n[full]  :  All terms are known: gate + C1,C2,Cm\n'

print(skewLabel)
print(fullLabel)

invalidChoice = True
while invalidChoice:
    choice = str(raw_input("[skew] or [full]? "))

    if choice == 'skew' or choice == 'full':
        invalidChoice = False
    else:
        print('Invalid entry. Try again.')

# temporary
# dotCapacitance = np.array([[49.1, -5.4],[-5.4, 54.7]])
# gateCapacitance = np.array([[-5.8, -0.6],[-0.6, -6.9]])


if choice == 'skew':
    # Get values from input CSV
    Cg1d1 = data['C_g1_d1'].values[0]
    Cg2d2 = data['C_g2_d2'].values[0]
    Cg1d2 = data['C_g1_d2'].values[0]
    Cg2d1 = data['C_g2_d1'].values[0]

    G = np.array([[1, Cg1d2/Cg2d2],[Cg2d1/Cg1d1, 1]])
    Gin = np.linalg.inv(G)

elif choice == 'full':
    # get values from csv
    Cg1d1 = data['C_g1_d1'].values[0]
    Cg2d2 = data['C_g2_d2'].values[0]
    Cg1d2 = data['C_g1_d2'].values[0]
    Cg2d1 = data['C_g2_d1'].values[0]
    C1 = data['C1'].values[0]
    C2 = data['C2'].values[0]
    Cm = data['Cm'].values[0]

    dotCapacitance = np.array([[C1, -1*Cm],[-1*Cm, C2]])
    gateCapacitance = -1*np.array([[Cg1d1, Cg2d1], [Cg1d2, Cg2d2]])

    G = np.dot(np.linalg.inv(dotCapacitance),gateCapacitance)
    # determinant instead?
    G = G/G[0][0]
    Gin = np.linalg.inv(G)

else:
    G = np.identity(2)
    Gin = G

# Get sweep input
invalidSweep = True
while invalidSweep:
    print('\nEnter intended sweep range for first gate voltage:')
    initial1 = raw_input("Initial: ")
    final1 = raw_input("Final: ")

    print('\nEnter intended sweep range for second gate voltage:')
    initial2 = raw_input("Initial: ")
    final2 = raw_input("Final: ")

    invalidSweep = False
    try:
        initial1 = float(initial1)
        initial2 = float(initial2)
        final1 = float(final1)
        final2 = float(final2)
    except ValueError as e:
        print('\nSweep range input not valid input.\n')
        invalidSweep = True

# Get virtual sweep range to match real sweep
virtInitial = np.dot(G,[initial1,initial2])
virtFinal = np.dot(G,[final1,final2])

print('\nReal gate voltage to virtual gate matrix:')
print(G)
print('\nVirtual gate to real gate voltage matrix:')
print(Gin)

# Formatting for printing equations
eqs = ['']*len(Gin)

virtual_variable_naming = 'VirtVar'
# virtual_variable_naming = 'u'
virtual_variable_index = 1
# virtual_variable_index = 0

for i,row in enumerate(Gin):
    for j,val in enumerate(row):
        eqs[i] = eqs[i] + '({0})*{1}{2}+'.format(val,virtual_variable_naming,j+virtual_variable_index)
    # cut off trailing '+'
    eqs[i] = eqs[i][:-1]


message  = '\n*******************************************************************************\n'
message += 'Virtual Setup instructions:\n'
message += 'To have the virtual gate sweeps cover the same operating region of the charge\n'
message += 'stability diagram, use the following parameters in the \'Virtual Setup\' application:\n\n'

message += 'Virtual Gate 1:\n'
message += '	name: vg1\n'
message += '	Initial: {0}\n'.format(virtInitial[0])
message += '	Final: {0}\n'.format(virtFinal[0])
message += '	Order: 1 \n\n'

message += 'Virtual Gate 2:\n'
message += '	name: vg2\n'
message += '	Initial: {0}\n'.format(virtInitial[1])
message += '	Final: {0}\n'.format(virtFinal[1])
message += '	Order: 2 \n\n'

message += 'Dependent variable expressions:\n\n'
message += '	gate1 = {0} \n\n'.format(eqs[0])
message += '	gate2 = {0} \n\n'.format(eqs[1])
message += '*******************************************************************************'

print(message)
