# -*- coding: utf-8 -*-
# Coloring Graph using HillClimbing algorithm and Best-First-Search algorithms
# Author: Luan Pham
# Date: 00:53 01-May-18
import copy
import random
import queue as Q

# GRAPH IS PRESENTED AS A ADJACENT MATRIX
ADJMAX =[
        [0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0],
        [0, 0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0],
        [0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1],
        [0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 0],
        [0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 1],
        [0, 1, 0, 1, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1],
        [1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 0, 0, 1, 1, 1, 0, 1],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 0, 1],
        [1, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1],
        [1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 0, 1, 0],
        [1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 0],
        [0, 1, 1, 0, 1, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0, 0, 1],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1],
        [1, 0, 0, 1, 0, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1, 0, 0, 1],
        [0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0],
        [0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 0, 0, 1, 1, 1, 1, 0, 1, 0, 0, 1],
        [0, 1, 1, 1, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0],
        [0, 1, 0, 0, 0, 1, 0, 1, 1, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 1, 0],
        [1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 0, 1, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0],
        [1, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0],
        [0, 1, 0, 1, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 1],
        [1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0],
        [0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1],
        [0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 1],
        [0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 0],
        [0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0, 1, 0],
        [0, 0, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 1, 1, 1],
        [1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 1, 1, 0, 1, 0, 0, 1],
        [1, 1, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0],
        [0, 0, 1, 0, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 0],
         ]	

# COLOR IS NUMBERED, INIT IS A ARRAY WITH SAME COLOR.
COLOR_ARRAY = [
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
        ]

# for test 
ADJMIN =[
        [0, 1, 1, 1],
        [1, 0, 1, 0],
        [1, 1, 0, 1],
        [1, 1, 1, 0],
        ]	

COLOR_TEST = [0, 0, 0, 0]

COLOR_SET = {
        0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20
        }

class myGraph: 
    def __init__(self, color, graph):   # DONE  INIT
        self.graph = graph
        self.color = color
        self.error = 999
    
    def displayColor(self):             # DONE  DISPLAY COLOR ARRAY
        print(self.color)

    def evaluate(self,color):           # DONE EVALUATE FOR HII CLIMBING
        evalError = 0
        for i in range(len(self.graph)):
            for j in range(len(self.graph)):
                if(self.graph[i][j] == 1 and color[i] == color[j]):
                    evalError += 1
        return evalError 

    def HillClimbing(self):             # DONE HILL CLIMBING
    
        # LOOP UNTIL SOLUTION IS FOUND
        while self.error != 0:

            # traverse color list
            for index in range(len(self.color)):
                    
                fillError = self.error
                fillColor = self.color[index]

                # traverse set color to find new operator
                for iColor in COLOR_SET:
                    
                    newColor = copy.copy(self.color)
                    
                    # SELECT AND APPLY A NEW OPERATOR
                    newColor[index] = iColor 
                    
                    # EVALUATE NEW OPERATOR
                    newError = self.evaluate(newColor)
                    if newError < fillError:
                        fillColor = iColor
                        fillError = newError
                        break
                
                
                # fill (Color,Error) is better right now
                self.color[index] = fillColor
                self.error       = fillError 
    
                # IF GOAL --> QUIT 
                if self.error == 0:
                    print("SUCCESS")
                    print("Graph is coloring by", len(set(self.color)), "color:")
                    self.displayColor()
                    break
    
    def BestFirst(self):
    
        # Create a priority Queue
        q = Q.PriorityQueue()
        
        # insert start node
        q.put((self.error, 0, self.color))
        
        # Until Priority Queue is empty
        while q.not_empty:

            item = q.get()   # item = (error, index, colorlist)
            iError = item[0]
            iIndex = item[1]
            iColor = item[2]

            # if Goal
            if iError == 0: 
           
                print("SUCCESS")
                print("Graph is coloring by", len(set(iColor)), "color:")
                print(iColor)
                break
            
            EvalMin = 999
            # put all new operator to queue
            for colorItem in set(iColor):
               
                newColor = copy.copy(iColor)
                newColor[iIndex] = colorItem 
                newError = self.evaluate(newColor)
                if newError < EvalMin:
                    EvalMin = newError

                q.put((newError, iIndex + 1, newColor))

            # try to insert one new
            newColor = copy.copy(iColor)
            newColor[iIndex] = max(iColor) + 1
            newError = self.evaluate(newColor)
            if newError < EvalMin:
                q.put((newError, iIndex + 1, newColor))
                    



# PROGRAM START HERE
print("Program start")
print("===================")
print("Coloring graph with BestFirst Search algorithms")
a = myGraph(COLOR_ARRAY, ADJMAX)
a.displayColor()
a.BestFirst()
print("===================")

print("Coloring graph with Hill Climbing algorithms")
b = myGraph(COLOR_ARRAY, ADJMAX)
b.displayColor()
b.HillClimbing()

















