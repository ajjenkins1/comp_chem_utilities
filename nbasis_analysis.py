import basis_set_exchange as bse
from ase.data import atomic_numbers, chemical_symbols

def calculate_nbasis(elements,basis_set):
    # Count the number of basis functions
    nbasis = sum(len(basis_set['elements'][str(el)]['electron_shells']) for el in elements)
    return nbasis

def is_integer(s):
    try:
        int(s)
        return True
    except ValueError:
        return False

# function to get the list of elements from an xyz file

def get_elements_from_xyz_file(geomfile):
    elements_numbers = []
    elements_symbols = []
    
    f=open(geomfile,"r")
    # first two lines of XYZ file are not the molecule
    lines = f.readlines()[2:]
    for x in lines:
        # get Atomic Symbol/Number 
        AtSym = x.split()[0] #('  ')[0]
        #In the XYZ file, Are they Symbols or Atomic Numbers?
        if is_integer(AtSym):
            atomic_number = int(AtSym)
            symbol = chemical_symbols[atomic_number]
            #print(f"{AtSym} is an atomic number. Symbol: {symbol}")
        else:
            symbol = AtSym
            atomic_number = atomic_numbers.get(symbol, None)
            #if atomic_number:
            #    #print(f"{AtSym} is an atomic symbol. Atomic number: {atomic_number}")
            #else:
            #    #print(f"{AtSym} is not a valid atomic symbol or number.")
        elements_numbers.append(atomic_number)
        elements_symbols.append(symbol)
    f.close
    return elements_numbers, elements_symbols

# Does this string contain the chars in list?
def ab(lelist, lestring):
    test = [e in lestring for e in lelist] # if e in lestring]
    if(test == [True]):
        testL = True
    else:
        testL = False
    return testL

# Does this string contain the chars in list?
def ab2(lelist, lestring):
    testL = False
    return (lelist in lestring)

# Find the highest angular momentum function for this basis+atom combo
def find_atom_highest_ang_mom(basis_set,atom):
    bs_str = bse.get_basis(basis_set, elements=atom, fmt='nwchem', header=False)
    l1 = bs_str.split('\n')[1]
    l2 = l1.split('[')
    if(ab(['h'],l2[1]) == True):
        high_ang_mom = 'h'
    elif(ab(['g'],l2[1]) == True):
        high_ang_mom = 'g'
    elif(ab(['f'],l2[1]) == True):
        high_ang_mom = 'f'
    elif(ab(['d'],l2[1]) == True):
        high_ang_mom = 'd'
    elif(ab(['p'],l2[1]) == True):
        high_ang_mom = 'p'
    elif(ab(['s'],l2[1]) == True):
        high_ang_mom = 's'
    else:
        high_ang_mom = 'X'
    #print(i," Highest ang mom: ", high_ang_mom)
    return high_ang_mom

# Does this atom+basis use an ECP?
def uses_ecp(basis_set,atoms):
    # just look if 'ECP' is printed in basis set definition
    bs_str = bse.get_basis(basis_set, elements=atoms, fmt='nwchem', header=False)
    return ab2("ECP",bs_str)

# Predefined mapping of atomic numbers to periods (rows)
periods = {
    1: 1, 2: 1,
    3: 2, 4: 2, 5: 2, 6: 2, 7: 2, 8: 2, 9: 2, 10: 2,
    11: 3, 12: 3, 13: 3, 14: 3, 15: 3, 16: 3, 17: 3, 18: 3,
    19: 4, 20: 4, 21: 4, 22: 4, 23: 4, 24: 4, 25: 4, 26: 4, 27: 4, 28: 4, 29: 4, 30: 4, 31: 4, 32: 4, 33: 4, 34: 4, 35: 4, 36: 4,
    37: 5, 38: 5, 39: 5, 40: 5, 41: 5, 42: 5, 43: 5, 44: 5, 45: 5, 46: 5, 47: 5, 48: 5, 49: 5, 50: 5, 51: 5, 52: 5, 53: 5, 54: 5,
    55: 6, 56: 6, 57: 6, 58: 6, 59: 6, 60: 6, 61: 6, 62: 6, 63: 6, 64: 6, 65: 6, 66: 6, 67: 6, 68: 6, 69: 6, 70: 6, 71: 6, 72: 6, 73: 6, 74: 6, 75: 6, 76: 6, 77: 6, 78: 6, 79: 6, 80: 6, 81: 6, 82: 6, 83: 6, 84: 6, 85: 6, 86: 6,
    87: 7, 88: 7, 89: 7, 90: 7, 91: 7, 92: 7, 93: 7, 94: 7, 95: 7, 96: 7, 97: 7, 98: 7, 99: 7, 100: 7, 101: 7, 102: 7, 103: 7, 104: 7, 105: 7, 106: 7, 107: 7, 108: 7, 109: 7, 110: 7, 111: 7, 112: 7, 113: 7, 114: 7, 115: 7, 116: 7, 117: 7, 118: 7
}

def get_row_number(atomic_number):
    row_number = periods.get(atomic_number, None)
    return row_number

# Predefined mapping of atomic numbers to groups
groups = {
    1: 1, 2: 18,
    3: 1, 4: 2, 5: 13, 6: 14, 7: 15, 8: 16, 9: 17, 10: 18,
    11: 1, 12: 2, 13: 13, 14: 14, 15: 15, 16: 16, 17: 17, 18: 18,
    19: 1, 20: 2, 21: 3, 22: 4, 23: 5, 24: 6, 25: 7, 26: 8, 27: 9, 28: 10, 29: 11, 30: 12, 31: 13, 32: 14, 33: 15, 34: 16, 35: 17, 36: 18,
    37: 1, 38: 2, 39: 3, 40: 4, 41: 5, 42: 6, 43: 7, 44: 8, 45: 9, 46: 10, 47: 11, 48: 12, 49: 13, 50: 14, 51: 15, 52: 16, 53: 17, 54: 18,
    55: 1, 56: 2, 57: 0, 58: 0, 59: 0, 60: 0, 61: 0, 62: 0, 63: 0, 64: 0, 65: 0, 66: 0, 67: 0, 68: 0, 69: 0, 70: 0, 71: 0, 72: 4, 73: 5, 74: 6, 75: 7, 76: 8, 77: 9, 78: 10, 79: 11, 80: 12, 81: 13, 82: 14, 83: 15, 84: 16, 85: 17, 86: 18,
    87: 1, 88: 2, 89: 0, 90: 0, 91: 0, 92: 0, 93: 0, 94: 0, 95: 0, 96: 0, 97: 0, 98: 0, 99: 0, 100: 0, 101: 0, 102: 0, 103: 0, 104: 4, 105: 5, 106: 6, 107: 7, 108: 8, 109: 9, 110: 10, 111: 11, 112: 12, 113: 13, 114: 14, 115: 15, 116: 16, 117: 17, 118: 18
}

def get_group_number(atomic_number):
    group_number = groups.get(atomic_number, None)
    return group_number

## Transition Metals
# by row
TM4 = list(range(21, 30+1))
TM5 = list(range(39, 48+1))
TM6 = list(range(72, 80+1))
TM7 = list(range(104, 108+1))
transition_metals = TM4 + TM5 + TM6 + TM7

def is_transition_metal(atomic_number):
    return atomic_number in transition_metals

def contains_transition_metal(elements):
    # could call is_transition_metal for each element but this is quicker
     return any(set(elements) & set(transition_metals))


## f-block
lanthanides = list(range(57,71+1))
actinides = list(range(89,103+1))

def is_lanthanide(atomic_number):
    #return 57 <= atomic_number <= 71
    # wrapper to allow single atomnic number query
    # use: is_actinide(100)
    return contains_lanthanide( list([atomic_number]) )

def contains_lanthanide(elements):
    # could call is_transition_metal for each element but this is quicker
     return any(set(elements) & set(lanthanides))


def is_actinide(atomic_number):
    #return 89 <= atomic_number <= 103
    # wrapper to allow single atomnic number query
    # use: is_actinide(100)
    return contains_actinide( list([atomic_number]) )

def contains_actinide(elements):
    # could call is_transition_metal for each element but this is quicker
     return any(set(elements) & set(actinides))

#The contains ones actually work if there is a list of one element long


