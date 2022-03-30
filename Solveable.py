
from abc import ABC, abstractmethod

class Solvable(ABC):
    @abstractmethod
    def tangentStiffness(self):
        pass
    
    @abstractmethod
    def updateGeometry(self):
        pass

