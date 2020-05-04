# Trash-Route-Optimization-Thesis-Project

Project Description
This repository includes the codes used to calculate the shortest distance route for a multiple Traveling Salesperson model.
The model is solved with MATLAB's intlinprog() function for interger linear programming. It has been applied to the Arizona 
State University Tempe campus Dual Litter containers waste/recycling collection route. For any number of travelers, their routes 
altogether will stop at each of the input locations starting from and ending at the first node input, the depot, while minimizing 
overall distance. The input data includes: node locations (coordinates), distance adjacency matrix, and user defined variables: 
number of travelers, min and max stops in each traveler's tour.

Additionally, the compactor code contains new linear constraints that account for the mid-tour stops at waste/recycling compactors 
that the ASU travelers must take to unload vehicles as they become full after visiting an input number of node loactions, H.


Organization
Data: Excel or MATLAB data files store the necessary input data for functions.

Functions: Regular or compactor myintlinprog() functions were created to be used depending on whether the route necessitates 
           compactor stops. These will solve the optimal routes given input data and plot them in MATLAB.

Size: The intlinprog() function will take longer to solve for larger problems. As a result, the ASU problem size was reduced
      by dividing the campus into zones with the K-means clustering function and routes found in smaller areas.
