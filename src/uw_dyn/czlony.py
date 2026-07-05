# -*- coding: utf-8 -*-
# uw_dyn: dynamika 3D ukladow wieloczlonowych
# autor: Marcin Szalajski (praca magisterska, 2016; pakiet od 2026)

"""Czlony (bryly sztywne) ukladu."""

import numpy as np


class Czlon:
    """Czlon (bryla sztywna): numer, masa i tensor bezwladnosci w ukladzie ciala."""
    
    def __init__(self, numer, masa, tensor_bez):
        self.i = numer 
        self.m = masa
        self.J = tensor_bez
    
    # macierz masowa    
    def M(self):
        return self.m*np.eye(3)
    
    
