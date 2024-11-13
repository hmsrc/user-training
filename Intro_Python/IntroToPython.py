
# coding: utf-8

# # This is a training supplement for HMS-RC Intro to Python Training
# # THIS SCRIPT FILE IS NOT DESIGNED TO BE EXECUTED.


# - The most recent default version can be found at [https://github.com/hmsrc/user-training](https://github.com/hmsrc/user-training)
# 
# - This notebook is configured to run with the virtual env found at `/n/groups/rc-training/python/env/trainingenv` on O2

# # Let's Learn Python!

# # Basic Syntax
# 
# - Statements are terminated by newlines (e.g. enter key)
# - Examples:

# In[ ]:


a = 1
print(a)


# In[ ]:


print("hello world")


# # Nuances: Quotations
# 
# - In general, double and single quotes are interchangable, both for printing and argument passing

# In[ ]:


print("hello world")
print('hello world')


# In[ ]:


str1 = 'abcdefg'
str2 = "abcdefg"

print(str1)
print(str2)


# In[ ]:


str1 == str2



# # Nuances: Escaped Characters
# 
# - Sometimes, you'll need to escape (backslash) certain special characters to get them to display correctly when printing. For example, if you want to preserve double quotations in your string:

# In[ ]:

print(""hello world"")

# In[ ]:


print("\"hello world\"")


# - A list of escape characters can be found in [Python documentation](https://docs.python.org/3/reference/lexical_analysis.html#literals) and elsewhere on the internet

# # Comments
# 
# - Comments are used to make code more legible. In python, they are denoted with the octothorpe.

# In[ ]:


# python will ignore this line
b = 2
print(b)
c = 3 #python will assign 'c = 3'
print(c)
d = 4 #python ignores everything after the pound, including: e = 5
print(d)
#print(e)


# - Multi-line comments can either be manually be broken down, or held in a docstring with triple quotes:

# In[ ]:


"""
this is a
multi-line comment
"""

print("Does my multi-line comment print out?")


# # Basic Features
# - Data Structures, Looping, etc.

# # Data Types
# 
# - Python is dynamically typed; there is no need to declare (e.g. int n = num). You can also reassign:

# In[ ]:


number = 1
print(number)


# In[ ]:


number = 2
print(number)


# In[ ]:


str3 = 'abc'
print(str3)


# In[ ]:


str3 = 'def'
print(str3)


# In[ ]:


str3 = 1
print(str3)


# # Data Types, an aside
# 
# - However, strings (and some other stuff like tuples) are not *mutable*; you cannot modify its content. When you reassign the variable as seen previously, you're actually generating a new object. You can't modify the original object:

# In[ ]:


s = "abc"
print(s[0])

# In[ ]:


#s[0] = "o"


# - The above example was shamelessly lifted off [Stack Overflow](http://stackoverflow.com/questions/8056130/immutable-vs-mutable-types-python)

# # Back to data types
# 
# - Like most languages, Python has higher level data structures. Of note are lists (arrays) and dictionaries. Relatedly, there are also tuples and sets. Each structure fills different niches.

# # Lists
# 
# - Lists are the most versatile; they are the effective equivalent of arrays in other languages. They are ordered, and you can fetch specific indices. You can mix and match types, have duplicates, nest lists, etc.

# To generate a list:

# In[ ]:


lst = [] #initialize a list of indeterminate size
lst.append('a')
print(lst)


# In[ ]:


lst.extend(['b', 'c'])
print(lst)


# In[ ]:


lst2 = [None]*5 #initialize a list of predetermined size (an array)
print(lst2)


# In[ ]:


lst2[1] = 'foo'
print(lst2[1])
print(lst2)


# # Lists, cont.

# In[ ]:


lst.insert(2, 'new_element')
print(lst)


# In[ ]:


lst.pop(2)
print(lst)


# In[ ]:


lst.pop()
print(lst)


# # Basic List Manipulations
# 
# - lists are zero-indexed (count from 0), and you can reference positions like so:
# 

# In[ ]:


lst = [1, 2, 3, 4, 5]
print(lst[1])


# In[ ]:


print(lst[1:])


# In[ ]:


print(lst[1:3])


# In[ ]:


print(lst[:3])


# In[ ]:


print(lst[1::2])


# In[ ]:


print(lst[-1])


# In[ ]:


lst2 = [1, 2, 3, [4, 5]]
print(lst2[3][1])


# # Sets
# 
# - Sets are lists that do not allow duplicates. Sets are useful if you need to take unions, intersections, etc. 

# To generate sets:

# In[ ]:


st = set(lst)
print(st)


# In[ ]:


lst3 = [1, 1, 2, 3, 4, 4]
st2 = set(lst3)
print(st2)


# In[ ]:


st3 = {1, 2, 3}
print(st3)


# # Sample Set Operations

# - Here are some elementary operations you can perform with sets:

# In[ ]:


set1 = {1, 2, 3}
set2 = {3, 4, 5}
print(set1 | set2)


# In[ ]:


print(set1 & set2)


# In[ ]:


print(set1 - set2)


# In[ ]:


print(set2 - set1)


# In[ ]:


print(set1 ^ set2)


# In[ ]:


print(set2 ^ set1)


# - More information can be found in the [documentation](https://docs.python.org/3/library/stdtypes.html?highlight=sets#set-types-set-frozenset).

# ## Tuples
# 
# Tuples are lists, but immutable. Once created, they cannot be modified. If you know the size of your data structure, tuples are preferred over lists because Python will know exactly how much memory to allocate.

# To create tuples:

# In[ ]:


empty_tuple=()
print(empty_tuple)


# In[ ]:


one_element = 1,
print(one_element)


# In[ ]:


mixed_tuple = (1, 'a', [1, 'a'])
print(mixed_tuple)


# In[ ]:


lst = [1, 2, 3]
tuple_from_list = tuple(lst)
print(tuple_from_list)


# ## Dictionaries
# 
# Dictionaries are associative collections (maps). They are composed of key:value pairs (most of the time). 
# 

# To create a dictionary:

# In[ ]:


# careful not to confuse this with empty sets;
# to initialize an empty set, use set())
empty_dict = {}  
print(empty_dict)


# In[ ]:


dict1 = {1: 'a', 2: 'b', 3: 'c'}
print(dict1)


# In[ ]:


dict2 = dict([(1, 'a'), (2, 'b'), (3, 'c')])
print(dict2)


# In[ ]:


# only works if keys are strings
dict3 = dict(a=1, b=2, c=3) 
print(dict3)


# In[ ]:


# see that this doesn't work 
# because the keys are not strings
dict4 = dict(1=a, 2=b, 3=c)
# will get an error


# ## Dictionaries, continued
# 

# In[ ]:


# to refresh your memory, this what is in dict1
print(dict1)


# In[ ]:


# add to dictionary:
dict1[4] ='d'
print(dict1)


# In[ ]:


# add to dictionary, but overwrite as key exists
dict1[4] = 'e'
print(dict1)


# In[ ]:


# add multiple key value pairs to dictionary
dict1.update({5:'f', 6:'g'}) 
print(dict1)


# In[ ]:


# remove a key value pair from dictionary
dict1.pop(1)
print(dict1)


# In[ ]:


# can also use del to remove key value pair
del dict1[2]
# will give an error if you have already deleted it


# ## Basic dictionary manipulations
# 
# Dictionaries don't have indices to reference, but they do have keys and values.

# To interact with dictionaries:

# In[ ]:


dict1 = {1: 'a', 2: 'b', 3: 'c'}
print(dict1[1])


# In[ ]:


print(dict1.keys())


# In[ ]:


print(dict1.values())


# In[ ]:


dict1.get(4, "I am a default value")


# There's a lot more you can do with these data structures. If you find yourself wondering if x data structure can do y thing, feel free to search for an answer. Often, there will be a solution!

# ## Flow Control
# 
# for, while, if/else

# ## For Loops
# 
# For loops are very straightforward to implement.

# To iterate over a list:

# In[ ]:


things = ['apple', 'banana', 'cherry']
for thing in things:
    print(thing)


# ## An aside: scoping
# 
# Code within indented blocks is known to be restricted in _scope_ to that block. That is, anything that happens within that block of code does not necessarily persist outside that block of code. We'll revisit this momentarily...

# ## For Loops, continued
# 
# This will go through every item in your list (or tuple or whatever) in order. If you also wanted to fetch indices, you'd use:

# In[ ]:


for idx, thing in enumerate(things):
    print(idx, thing)


# ## For Loops, continued
# 
# - The `enumerate()` function is actually an _iterator_ that spits out a tuple consisting of (index, value) and assigns each to idx, name. This is called a **named tuple**.
# - The `for a in b` structure can be replaced with any generic construct (within reason). For a generic loop, you can do something like:

# In[ ]:


for i in range(5):
    print(i)


# ## An application of for loops: File I/O
# 
# - File input/output (I/O) in python is conventionally done via the creation of a _context manager_. You can think of this (i.e. the file handle and its associated operations) as a block of code that executes for as long as you need that file to stay open.
# 
# ______
# 
# 
# ```
# with open('file.txt', 'r') as f:
#     for line in f:
#         print(line)
# ```
# 
# ______
# 
# 
# - The code runs in the context of the open file, and automatically closes the file when done.
# - This is an example of scoping. We'll have another example when scripting.
# 

# ## File I/O continued (an aside)
# 
# The first argument of `open()` is the filename. The second one (the mode) determines how the file gets opened.
# 
# - Read the file, pass in "r". Read and write the file, pass in "r+"
# - Overwrite the file, pass in "w". Writing and reading the file, pass "w+" Overwrites the existing file if it already exists. 
# - If you want to append to the file, pass in "a"
# 
# ______
# 
# 
# ```
# with open('file.txt', 'r', encoding="utf-8"): as f:
#     f.readlines(2) # Read in N # lines
#     for line in f: # Do something with each line
#         x = f.readline()
# ```
# ______
# 
# ______
# 
# 
# ```
# import json, os
# with open('file.txt', 'w', encoding="utf-8"): as f:
#     for x in range(100):
#         f.writeline(x)
#     # write all env vars to file
#     f.write(json.dumps(os.environ, indent=2))
# ```
# ______
# 
# 
# 
# 
# If you are getting unexpected input when performing file IO, try adding the following:
# 
# ______
# 
# ```
# with open('file.txt', 'r', encoding="utf-8"): as f:​
# ```
# ______
# 

# ## While Loops
# 
# Similar in function to for loops, while loops are used to iterate, and are useful if you don't know how long to iterate for (e.g. searching for convergence, etc.). A generic while loop looks like this:

# In[ ]:


toggle = True

while toggle == True:
    toggle = False
    print(toggle)


# A very simplistic example, but the above while loop runs one iteration, then exits because `toggle` is no longer `True`.

# ## If/elif/else
# 
# `if/else` statements are useful when you need to handle different cases in your workflow. A generic if/elif/else statement structure:
# 
# _____
# 
# ```
# # take some input (let's say, a number)​
# if input == 0:       #if input is zero​
#     print('zero')​
# elif input % 2 == 0: #if input is otherwise even​
#     print('even')​
# else:                #if input is odd​
#     print ('odd')
# ```
# _____
# 
# You may use as many elifs as you require to solve your problem.

# ## An aside: try/except
# 
# Python has the interesting distinction of relatively lightweight exception handling. Exceptions are only computationally expensive if they trigger. For more information, [this page has a good primer on try/except](https://www.askpython.com/python/python-exception-handling) and [the Python documentation on try/except is located here.](https://docs.python.org/3/tutorial/errors.html)

# ## A few pertinent modules    
#     
# NumPy is the go-to module if you need to perform complex arithmetic or do matrix operations (manipulate data frames). It is much more efficient than using system Python utilities.
# 
# SciPy has a bunch of interesting functions that automate various aspects of scientific computing, like statistics and higher-level mathematics. Most of these are callable in one line. ([documentation](http://docs.scipy.org/doc/), includes NumPy)
# 
# Matplotlib is a basic plotting module. You can generate plots in real time (X11) or create them and write to files to view later. There is also a degree of customization afforded in plots. ([documentation](https://matplotlib.org/stable/users/index.html))
# 
# A super quick overview of some examples from each module follows.

# ### `scipy`

# In[ ]:


from scipy import constants        # `from scipy import *` will not behave as expected
print(constants.c)                 # speed of light
from scipy.stats import norm
# cumulative distribution function for normal continuous random variable
print(norm.cdf(5, 0, 3))           # P(x<5) if X~N(0,3)


# ### `numpy`

# In[ ]:


import numpy as np          # the general way; using `as` allows you to set an alias (if the name is too long to type)
cvalues = [25.3, 24.8, 26.9, 23.9]    # a typical python list of Celcius values
print(cvalues)
C = np.array(cvalues)       # create a numpy array
print(C)                    # note that they look identical
print(C * 9 / 5 + 32)       # convert to Fahrenheit using scalar multiplication
fvalues = [x*9/5 + 32 for x in cvalues]    # with a typical list, you need to loop over it and compute element-wise instead
print(fvalues)              # same result, less efficient/readable


# #### A brief look at arrays: NumPy data structures
# Note the formatting of the print below. This is indicative of the logic and illustrates why numpy is the preferred implementation of Python matrix operations. The SciPy documentation linked previously also includes NumPy documentation.

# In[ ]:


nested_list = [[3.4, 8.7, 9.9], [1.1, -7.8, -0.7], [4.1, 12.3, 4.8]]      # a python nested list
print(nested_list)
A = np.array ([ [3.4, 8.7, 9.9], [1.1, -7.8, -0.7], [4.1, 12.3, 4.8]])    # a numpy array (matrix)
print(A)


# ### `matplotlib`

# In[ ]:


import matplotlib.pyplot as plt


# In[ ]:


# the first line of code is a "magic command" from the iPython kernel:
# https://ipython.readthedocs.io/en/stable/interactive/magics.html#magic-matplotlib
# it "set up matplotlib to work interactively" and
# allows the plot to render below the code block in the notebook.
#get_ipython().run_line_magic('matplotlib', 'inline')
fig = plt.figure()
plt.plot([1,2,3,4])             # this will output some memory address
plt.ylabel('some numbers')  # another memory address
# in other contexts, you'll want to use plt.show() to see your plot 
plt.show()


# Let's go through a (quick) full example with some fake data. We'll demo a rudimentary K-means clustering implementation using the above packages. Remember, the focus here is not on what the functions do (though they'll be defined), but the power that these packages hold together. If you have questions about the imported modules and functions leveraged, we recommend consulting the documentation (linked previously).
# 
# The process for following example has been co-opted from https://www.tutorialspoint.com/scipy/scipy_cluster.htm, so if you'd like a bit more information about this specific code stack, feel free to look over there.

# In[ ]:


# first, let's generate some data.
from numpy import vstack,array
from numpy.random import rand

# data generation with three features
data = vstack((rand(100,3) + array([.5,.5,.5]),rand(100,3)))
print(data.shape) # gives "shape" of the numpy array
                  # as this is a 2D array, it is the num rows x num columns
print(data)


# In[ ]:


# whiten the data
# https://en.wikipedia.org/wiki/Whitening_transformation
from scipy.cluster.vq import whiten

data = whiten(data)
print(data)


# In[ ]:


#get_ipython().run_line_magic('matplotlib', 'ipympl')

# visualize the data:
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure(1)
ax = fig.add_subplot(111, projection='3d')
ax.scatter(data[:,0], data[:,1], data[:,2])
plt.show()


# In[ ]:


# perform k-means with 3 clusters
from scipy.cluster.vq import kmeans

#_ is a placeholder. kmeans returns 2 objects
# we don't care about the second returned value here
centroids,_ = kmeans(data,3)
print(centroids)


# In[ ]:


#get_ipython().run_line_magic('matplotlib', 'ipympl')

# superimpose centroids/recreate plot
fig2 = plt.figure(2)
ax2 = fig2.add_subplot(111, projection='3d')
ax2.scatter(data[:,0], data[:,1], data[:,2])
ax2.scatter(centroids[:,0], centroids[:,1], centroids[:,2], s=300)
plt.show()

# In[ ]:


# distribute data to each cluster
from scipy.cluster.vq import vq

# using _ again to discard returned results we don't need
assignments,_ = vq(data, centroids)
print(assignments)


# In[ ]:


#get_ipython().run_line_magic('matplotlib', 'ipympl')

# regenerate the plot with clusters color-coded
from matplotlib.colors import ListedColormap

colors = ['red', 'green', 'blue'] # colors of 3 clusters

fig3 = plt.figure(3)
ax3 = fig3.add_subplot(111, projection='3d')
ax3.scatter(data[:,0], data[:,1], data[:,2], c=assignments, cmap=ListedColormap(colors))
ax3.scatter(centroids[:,0], centroids[:,1], centroids[:,2], marker='*', s=300, c=colors, cmap=ListedColormap(colors))
plt.show()


# # Up Next: Object Oriented Programming
# 
# - Please switch back to the powerpoint.
