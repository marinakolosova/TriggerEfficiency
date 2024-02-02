import math
import itertools
import logging
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

def deltaPhi(phi1, phi2):
    try:
        dphi = phi1 - phi2
    except TypeError:
        dphi = phi1.phi - phi2.phi
    while dphi > math.pi:
        dphi -= 2 * math.pi
    while dphi < -math.pi:
        dphi += 2 * math.pi
    return dphi

def deltaR2(eta1, phi1, eta2=None, phi2=None):
    if eta2 is None:
        a, b = eta1, phi1
        return deltaR2(a.eta, a.phi, b.eta, b.phi)
    else:
        deta = eta1 - eta2
        dphi = deltaPhi(phi1, phi2)
        return deta * deta + dphi * dphi

def deltaR(eta1, phi1, eta2=None, phi2=None):
    return math.sqrt(deltaR2(eta1, phi1, eta2, phi2))

def deltaEta(obj1, obj2):
    return abs(obj1.eta - obj2.eta)
