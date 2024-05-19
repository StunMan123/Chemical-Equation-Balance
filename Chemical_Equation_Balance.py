import re
import sympy
import scipy.optimize
import numpy as np
from sympy import Add
from sympy import Matrix, nsimplify, Rational
from scipy.optimize import least_squares

def addToMatrix(element, index, count, side):
    global elementList, elementMatrix
    if index == len(elementMatrix):
        elementMatrix.append([0] * len(elementList))
    if element not in elementList:
        elementList.append(element)
        for row in elementMatrix:
            row.append(0)
    column = elementList.index(element)
    elementMatrix[index][column] += count * side


def findElements(segment, index, multiplier, side):
    global elementList, elementMatrix
    elementsAndNumbers = re.split('([A-Z][a-z]?)', segment)
    i = 0
    while i < len(elementsAndNumbers) - 1:  # last element always blank
        i += 1
        if len(elementsAndNumbers[i]) > 0:
            if i + 1 < len(elementsAndNumbers) and elementsAndNumbers[i + 1].isdigit():
                count = int(elementsAndNumbers[i + 1]) * multiplier
                addToMatrix(elementsAndNumbers[i], index, count, side)
                i += 1
            else:
                addToMatrix(elementsAndNumbers[i], index, multiplier, side)


def compoundDecipher(compound, index, side):
    global elementList, elementMatrix
    segments = re.split('(\([A-Za-z0-9]*\)[0-9]*)', compound)
    for segment in segments:
        if segment.startswith("("):
            segment = re.split('\)([0-9]*)', segment)
            multiplier = int(segment[1])
            segment = segment[0][1:]
        else:
            multiplier = 1
        findElements(segment, index, multiplier, side)

def balance_equation(reactants, products):
    global elementList, elementMatrix
    # Reset the global variables to start fresh each time the function is called
    elementList = []
    elementMatrix = []

    # Fill the element matrix for reactants and products
    for i, reactant in enumerate(reactants):
        compoundDecipher(reactant, i, 1)
    for i, product in enumerate(products):
        compoundDecipher(product, i + len(reactants), -1)

    # Ensure elementMatrix is a 2D numpy array
    elementMatrix = np.array(elementMatrix, dtype=float)

    # Transpose the matrix to align rows with elements and columns with compounds
    elementMatrix = elementMatrix.T

    # Define the objective function to calculate the residuals
    def objective_function(coeffs):
        # Calculate the product of coefficients and the element matrix
        return np.dot(elementMatrix, coeffs)

    # Set up the initial guess for the coefficients as ones
    initial_guess = np.ones(len(reactants) + len(products))

    # Use least squares optimization to find the best coefficients
    result = least_squares(objective_function, initial_guess, bounds=(0, np.inf))

    if not result.success:
        raise ValueError("Optimization failed: " + result.message)

    # Scale coefficients to make the first one equal to 1
    coefficients = result.x / result.x[0]

    return coefficients


# Read reactants and products from user input
print("Please input your reactants, this is case sensitive")
reactants = input("Reactants: ")
print("Please input your products, this is case sensitive")
products = input("Products: ")

reactants = reactants.replace(' ', '').split("+")
products = products.replace(' ', '').split("+")

elementList = []
elementMatrix = []

# Solve the equation
coefficients = balance_equation(reactants, products)

# Generate the output
output = ""
for i, reactant in enumerate(reactants):
    coeff = np.format_float_positional(coefficients[i], precision=4, unique=False, fractional=False, trim='k')
    output += (coeff if coeff != '1' else '') + reactant
    if i < len(reactants) - 1:
        output += " + "
output += " -> "
for i, product in enumerate(products):
    coeff = np.format_float_positional(coefficients[i + len(reactants)], precision=4, unique=False, fractional=False, trim='k')
    output += (coeff if coeff != '1' else '') + product
    if i < len(products) - 1:
        output += " + "

print(output)


#C6H10O5+C5H8O4+C18H14O7+H2
#C12H26+C6H6+C9H12+C7H8+H2O

#C8H18+C9H20+C10H22+C11H24+C12H26+C13H28+C14H30+C15H32+C16H34+H2O

#CH4+C2H6+C3H8+C4H10+C5H12+C6H14+C7H16+C8H18+C9H20+C10H22+C11H24+C12H26+C13H28+C14H30+C15H32+C16H34+H2O
#C6H6O3+C5H4O2+C5H8O2
#C18H14O7(lignin)
#C6H10O5(cellulose)
#C5H8O4(hemicelllulose)
#C6H6O3(5-HMF)
#C10H10O4(5-HMF dimer)
#C14H14O5(5-HMF trimer)
#C18H18O6(5-HMF tetramer)
#C4H4O(furan)
#C5H4O2(furfural)
#C9H8O3(furfural dimer)
#C13H12O4(furfural trimer)

#C6H6(benzene)
#C9H12(cumene)
#C7H8(toluene)

#C6H12O6+NH4NO3
#C4H8N2O3+CH1.8O0.5N0.2+CO2+H2O