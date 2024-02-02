#!/usr/bin/env python3
import os, sys
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
from importlib import import_module

import array
from helpers.utils import deltaPhi, deltaR

#importing tools from nanoAOD processing set up to store the ratio histograms in a root file
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

class TrigHadAnalysis(Module):
    def __init__(self):
        self.writeHistFile=True
        self.reference_paths=reference_paths
        self.signal_paths=signal_paths
        
    def beginJob(self,histFile=None,histDirName=None):
        Module.beginJob(self,histFile,histDirName)

        self.h_passreftrig = ROOT.TH1F("h_passreftrig" , "; passed ref trigger", 2, 0. , 2.)
        self.h_dR_AK8_mu = ROOT.TH1F("h_dR_AK8_mu", "; #Delta R(AK8, #mu);Efficiency", 50, 0, 5.0)

        self.labels = {}
        self.labels["pt"] = ';AK8 jet p_{T} [GeV];Efficiency'
        self.labels["mSD"]= ';AK8 jet m_{SD} [GeV];Efficiency'
        self.labels["eta"]= ';AK8 jet #eta;Effciency'
        
        self.binning = {}
        self.binning["pt"] = [200, 250, 300, 350, 400, 425, 450, 475, 500, 550, 600, 650, 700, 800, 900, 1000]
        self.binning["mSD"] = [0, 10, 20, 30, 35, 40, 50, 60, 70, 80, 90, 100, 125, 150, 175, 200, 225, 250, 300, 400]
        self.binning["eta"] = [-2.4, -2.2, -2.0, -1.8, -1.6, -1.4, -1.2, -1.0, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4]
        
        self.variables = ["pt", "mSD", "eta"]
        self.nanoVars = {"pt": "pt",
                         "mSD": "msoftdrop",
                         "eta": "eta",
                         "phi": "phi"}
        
        self.hList = {}

        self.hList['h_AK8_mSD_vs_pt']= ROOT.TH2F('h_AK8_mSD_vs_pt', ';p_{T} [GeV];m_{SD} [GeV];Efficiency', len(self.binning["pt"])-1,  array.array("d", self.binning["pt"]), len(self.binning["mSD"])-1, array.array("d", self.binning["mSD"]))
        self.hList['h_AK8_mSD_vs_pt_passTrgOR'] = ROOT.TH2F('h_AK8_mSD_vs_pt_passTrgOR', ';p_{T} [GeV];m_{SD} [GeV];Efficiency', len(self.binning["pt"])-1, array.array("d", self.binning["pt"]), len(self.binning["mSD"])-1, array.array("d", self.binning["mSD"]))

        for var in self.variables:
            self.hList[f'h_AK8_{var}'] = ROOT.TH1F(f'h_AK8_{var}', self.labels[var], len(self.binning[var])-1, array.array("d", self.binning[var]))
            self.hList[f'h_AK8_{var}_plateauMSD'] = ROOT.TH1F(f'h_AK8_{var}_plateauMSD', self.labels[var], len(self.binning[var])-1, array.array("d", self.binning[var]))
            self.hList[f'h_AK8_{var}_plateauPt'] = ROOT.TH1F(f'h_AK8_{var}_plateauPt', self.labels[var], len(self.binning[var])-1, array.array("d", self.binning[var]))
            self.hList[f'h_AK8_{var}_plateauMSDPt'] = ROOT.TH1F(f'h_AK8_{var}_plateauMSDPt', self.labels[var], len(self.binning[var])-1, array.array("d", self.binning[var]))
            
            self.hList[f'h_AK8_{var}_passTrgOR'] = ROOT.TH1F(f'h_AK8_{var}_passTrgOR', self.labels[var], len(self.binning[var])-1, array.array("d", self.binning[var]))
            self.hList[f'h_AK8_{var}_passTrgOR_plateauMSD'] = ROOT.TH1F(f'h_AK8_{var}_passTrgOR_plateauMSD', self.labels[var], len(self.binning[var])-1, array.array("d", self.binning[var]))
            self.hList[f'h_AK8_{var}_passTrgOR_plateauPt'] = ROOT.TH1F(f'h_AK8_{var}_passTrgOR_plateauPt', self.labels[var], len(self.binning[var])-1, array.array("d", self.binning[var]))
            self.hList[f'h_AK8_{var}_passTrgOR_plateauMSDPt'] = ROOT.TH1F(f'h_AK8_{var}_passTrgOR_plateauMSDPt', self.labels[var], len(self.binning[var])-1, array.array("d", self.binning[var]))
            
            for path in self.signal_paths:
                self.hList[f'h_AK8_{var}_pass_HLT_{path}'] = ROOT.TH1F(f'h_AK8_{var}_pass_HLT_{path}', self.labels[var], len(self.binning[var])-1, array.array("d", self.binning[var]))
                self.hList[f'h_AK8_{var}_pass_HLT_{path}_plateauMSD']=ROOT.TH1F(f'h_AK8_{var}_pass_HLT_{path}_plateauMSD', self.labels[var], len(self.binning[var])-1, array.array("d", self.binning[var]))
                self.hList[f'h_AK8_{var}_pass_HLT_{path}_plateauPt']=ROOT.TH1F(f'h_AK8_{var}_pass_HLT_{path}_plateauPt', self.labels[var], len(self.binning[var])-1, array.array("d", self.binning[var]))
                self.hList[f'h_AK8_{var}_pass_HLT_{path}_plateauMSDPt']=ROOT.TH1F(f'h_AK8_{var}_pass_HLT_{path}_plateauMSDPt', self.labels[var], len(self.binning[var])-1, array.array("d", self.binning[var]))
                
        self.addObject(self.h_passreftrig)
        self.addObject(self.h_dR_AK8_mu)
        for h in self.hList:
            self.addObject(self.hList[h])

    def analyze(self, event):

        muons = Collection(event, "Muon")
        fatjets = Collection(event, "FatJet")
        hlt = Object(event, "HLT")
        
        # Check if event passes the reference trigger(s)
        refAccept=False
        for path in self.reference_paths:
            bit = getattr(hlt, path)
            if bit:
                refAccept = True

        # Save the bit of reference trigger and skim event
        self.h_passreftrig.Fill(refAccept)
        if not refAccept:
            return False

        # Add any offline selection here:
        selected_muon = []
        for mu in muons:
            if not mu.isGlobal:
                continue
            if mu.pt < 26:
                continue
            if abs(mu.eta) > 2.5:
                continue
            if abs(mu.dz) > 0.10:
                continue
            if abs(mu.dxy) > 0.05:
                continue
            if not abs(mu.looseId):
                continue
            selected_muon.append(mu)
        
        if len(selected_muon) == 0:
            return False

        #print("\n")
        #print("Selected muons: ", len(selected_muon))
        #for muIdx, muon in enumerate(selected_muon):
        #    print("muon ", muIdx, "   pt=", muon.pt, "   eta=", muon.eta, "   phi=", muon.phi)
        #print("Selected fatjets : ", len(fatjets))

        non_overlap_fatjets = []
        for fIdx, fatjet in enumerate(fatjets):
            # skip if fatjet overlaps with a muon
            if (fatjet.muonIdx3SJ != -1):
                dR = deltaR(fatjet.eta, fatjet.phi, muons[fatjet.muonIdx3SJ].eta, muons[fatjet.muonIdx3SJ].phi)
                self.h_dR_AK8_mu.Fill(dR)
            else:
                non_overlap_fatjets.append(fatjet)

        if len(non_overlap_fatjets) == 0:
            return False
        sorted_fatjets = sorted(non_overlap_fatjets, key=lambda x: x.pt, reverse=True)

        # Fill numerator
        for var in self.variables:
            self.hList[f'h_AK8_{var}'].Fill(getattr(sorted_fatjets[0], self.nanoVars[var]))
        self.hList['h_AK8_mSD_vs_pt'].Fill(sorted_fatjets[0].pt, sorted_fatjets[0].msoftdrop)

        signalOR = False
        for path in self.signal_paths:
            if getattr(hlt, path) == 1:
                signalOR = True
                for var in self.variables:
                    self.hList[f'h_AK8_{var}_pass_HLT_{path}'].Fill(getattr(sorted_fatjets[0], self.nanoVars[var]))
        
        if signalOR:
            self.hList['h_AK8_mSD_vs_pt_passTrgOR'].Fill(sorted_fatjets[0].pt, sorted_fatjets[0].msoftdrop)
            for var in self.variables:
                self.hList[f'h_AK8_{var}_passTrgOR'].Fill(getattr(sorted_fatjets[0], self.nanoVars[var]))

        # Require the AK8 jet to be on the plateau of the soft-drop mass:
        if sorted_fatjets[0].msoftdrop > 50:
            for var in self.variables:
                self.hList[f'h_AK8_{var}_plateauMSD'].Fill(getattr(sorted_fatjets[0], self.nanoVars[var]))
                if signalOR:
                    self.hList[f'h_AK8_{var}_passTrgOR_plateauMSD'].Fill(getattr(sorted_fatjets[0], self.nanoVars[var]))
                for path in self.signal_paths:
                    if getattr(hlt, path) == 1:
                        self.hList[f'h_AK8_{var}_pass_HLT_{path}_plateauMSD'].Fill(getattr(sorted_fatjets[0], self.nanoVars[var]))

        # Requre the AK8 jet to be on the plateau of the pT:
        if sorted_fatjets[0].pt > 480:
            for var in self.variables:
                self.hList[f'h_AK8_{var}_plateauPt'].Fill(getattr(sorted_fatjets[0], self.nanoVars[var]))
                if signalOR:
                    self.hList[f'h_AK8_{var}_passTrgOR_plateauPt'].Fill(getattr(sorted_fatjets[0], self.nanoVars[var]))
                for path in self.signal_paths:
                    if getattr(hlt, path) == 1:
                        self.hList[f'h_AK8_{var}_pass_HLT_{path}_plateauPt'].Fill(getattr(sorted_fatjets[0], self.nanoVars[var]))
                
        # Require the AK8 jet to be on the plateau of both legs:
        if sorted_fatjets[0].pt > 480 and sorted_fatjets[0].msoftdrop > 50:
            for var in self.variables:
                self.hList[f'h_AK8_{var}_plateauMSDPt'].Fill(getattr(sorted_fatjets[0], self.nanoVars[var]))
                if signalOR:
                    self.hList[f'h_AK8_{var}_passTrgOR_plateauMSDPt'].Fill(getattr(sorted_fatjets[0], self.nanoVars[var]))
                for path in self.signal_paths:
                    if getattr(hlt, path) == 1:
                        self.hList[f'h_AK8_{var}_pass_HLT_{path}_plateauMSDPt'].Fill(getattr(sorted_fatjets[0], self.nanoVars[var]))

        return True

preselection="(HLT_Mu50 == 1 || HLT_IsoMu24 == 1) && (Sum$(FatJet_pt > 200 && abs(FatJet_eta)<2.5) > 0)"
files=[
    "/eos/cms/store/group/dpg_trigger/comm_trigger/TriggerStudiesGroup/STEAM/hlt_tutorial/Muon2023C/03db0efb-09c4-4f94-b341-52e2e0947da5.root",
    "/eos/cms/store/group/dpg_trigger/comm_trigger/TriggerStudiesGroup/STEAM/hlt_tutorial/Muon2023C/0b30c7f6-0508-4a5e-9018-c9c95cd0b3bf.root",
    "/eos/cms/store/group/dpg_trigger/comm_trigger/TriggerStudiesGroup/STEAM/hlt_tutorial/Muon2023C/13390fda-f51d-49f1-853e-f1e6d40fba91.root",
    "/eos/cms/store/group/dpg_trigger/comm_trigger/TriggerStudiesGroup/STEAM/hlt_tutorial/Muon2023C/13fe7693-b3a8-4c1e-b378-7df0d6d01e77.root",
    "/eos/cms/store/group/dpg_trigger/comm_trigger/TriggerStudiesGroup/STEAM/hlt_tutorial/Muon2023C/34568e2a-a4d2-4a5f-8e7a-bbc15d842894.root",
    "/eos/cms/store/group/dpg_trigger/comm_trigger/TriggerStudiesGroup/STEAM/hlt_tutorial/Muon2023C/45fa89ac-f3f4-493a-bd77-fff590e507b1.root",
    "/eos/cms/store/group/dpg_trigger/comm_trigger/TriggerStudiesGroup/STEAM/hlt_tutorial/Muon2023C/48dbe055-58b3-40e4-9970-71fff20a8438.root",
    "/eos/cms/store/group/dpg_trigger/comm_trigger/TriggerStudiesGroup/STEAM/hlt_tutorial/Muon2023C/50eb4ba0-6af8-4563-a001-688ed0ae919a.root",
    "/eos/cms/store/group/dpg_trigger/comm_trigger/TriggerStudiesGroup/STEAM/hlt_tutorial/Muon2023C/5303c703-9856-4da2-8972-7a8cbce039b7.root",
    "/eos/cms/store/group/dpg_trigger/comm_trigger/TriggerStudiesGroup/STEAM/hlt_tutorial/Muon2023C/58b6f331-8fdc-400f-adfb-14c79a40567b.root",
    "/eos/cms/store/group/dpg_trigger/comm_trigger/TriggerStudiesGroup/STEAM/hlt_tutorial/Muon2023C/5b440dc4-2037-42a0-a97e-e58b64db0674.root",
    "/eos/cms/store/group/dpg_trigger/comm_trigger/TriggerStudiesGroup/STEAM/hlt_tutorial/Muon2023C/6c704034-88a5-4fc7-804d-1f2a3c88257e.root",
    "/eos/cms/store/group/dpg_trigger/comm_trigger/TriggerStudiesGroup/STEAM/hlt_tutorial/Muon2023C/808f4b4c-a7b9-4d13-ae1c-60575a5a1df8.root",
    "/eos/cms/store/group/dpg_trigger/comm_trigger/TriggerStudiesGroup/STEAM/hlt_tutorial/Muon2023C/8091b148-c47c-4634-af07-611e5f42e0a2.root",
    "/eos/cms/store/group/dpg_trigger/comm_trigger/TriggerStudiesGroup/STEAM/hlt_tutorial/Muon2023C/827c9945-8689-4598-a7b2-f94242e30498.root",
    "/eos/cms/store/group/dpg_trigger/comm_trigger/TriggerStudiesGroup/STEAM/hlt_tutorial/Muon2023C/8acd3808-bf6c-4d82-b0a2-4271efa9e274.root",
    "/eos/cms/store/group/dpg_trigger/comm_trigger/TriggerStudiesGroup/STEAM/hlt_tutorial/Muon2023C/8ec6ea95-591d-489f-92d3-4e6d9ceef5d9.root",
    "/eos/cms/store/group/dpg_trigger/comm_trigger/TriggerStudiesGroup/STEAM/hlt_tutorial/Muon2023C/9b57ac11-dbfa-414d-9dfe-66bf1d4dad1b.root",
    "/eos/cms/store/group/dpg_trigger/comm_trigger/TriggerStudiesGroup/STEAM/hlt_tutorial/Muon2023C/a2cbc9f4-8fa6-46fd-8e39-dedafecefeea.root",
    "/eos/cms/store/group/dpg_trigger/comm_trigger/TriggerStudiesGroup/STEAM/hlt_tutorial/Muon2023C/ada87593-bfc3-4f72-9272-22a05b0a582c.root",
    "/eos/cms/store/group/dpg_trigger/comm_trigger/TriggerStudiesGroup/STEAM/hlt_tutorial/Muon2023C/de7ca17b-e7cf-4fa5-81d4-e4ef3bc4a3de.root",
    "/eos/cms/store/group/dpg_trigger/comm_trigger/TriggerStudiesGroup/STEAM/hlt_tutorial/Muon2023C/bee92b12-e39d-4490-a641-b13655a5ed4f.root",
    "/eos/cms/store/group/dpg_trigger/comm_trigger/TriggerStudiesGroup/STEAM/hlt_tutorial/Muon2023C/e4d3e2e4-208f-41ad-ad1a-be4981457604.root",
    "/eos/cms/store/group/dpg_trigger/comm_trigger/TriggerStudiesGroup/STEAM/hlt_tutorial/Muon2023C/f7fea811-fe25-445e-b326-eae525a2dfd6.root",
    "/eos/cms/store/group/dpg_trigger/comm_trigger/TriggerStudiesGroup/STEAM/hlt_tutorial/Muon2023C/01adecde-d165-4569-a650-3cccc0971c4f.root",
    "/eos/cms/store/group/dpg_trigger/comm_trigger/TriggerStudiesGroup/STEAM/hlt_tutorial/Muon2023C/04ea2edc-f5ee-48a6-8420-27d57e2bc520.root",
    "/eos/cms/store/group/dpg_trigger/comm_trigger/TriggerStudiesGroup/STEAM/hlt_tutorial/Muon2023C/0d36eba2-dd59-4af5-8fb2-f883eefba4d4.root",
    "/eos/cms/store/group/dpg_trigger/comm_trigger/TriggerStudiesGroup/STEAM/hlt_tutorial/Muon2023C/33b52c18-f88c-4b1a-a85e-5bfc372cc75f.root",
    "/eos/cms/store/group/dpg_trigger/comm_trigger/TriggerStudiesGroup/STEAM/hlt_tutorial/Muon2023C/4a9ef871-d7a3-4267-98fc-79445aeb0ac2.root",
    "/eos/cms/store/group/dpg_trigger/comm_trigger/TriggerStudiesGroup/STEAM/hlt_tutorial/Muon2023C/4dda7e4c-f811-4526-be3f-f7ab5901d23c.root",
    "/eos/cms/store/group/dpg_trigger/comm_trigger/TriggerStudiesGroup/STEAM/hlt_tutorial/Muon2023C/66421788-2e1b-4f5d-9ef9-d2ca91597090.root",
    "/eos/cms/store/group/dpg_trigger/comm_trigger/TriggerStudiesGroup/STEAM/hlt_tutorial/Muon2023C/68fcb303-070e-4416-90bd-0ad5d5aab7fa.root",
    "/eos/cms/store/group/dpg_trigger/comm_trigger/TriggerStudiesGroup/STEAM/hlt_tutorial/Muon2023C/6c1267bb-51ff-4844-b28e-325c578d3dad.root",
    "/eos/cms/store/group/dpg_trigger/comm_trigger/TriggerStudiesGroup/STEAM/hlt_tutorial/Muon2023C/750c89cc-9469-468f-9e94-4cb8648488dc.root",
    "/eos/cms/store/group/dpg_trigger/comm_trigger/TriggerStudiesGroup/STEAM/hlt_tutorial/Muon2023C/916ef9fb-1400-4bc6-bd4b-368ded630492.root",
    "/eos/cms/store/group/dpg_trigger/comm_trigger/TriggerStudiesGroup/STEAM/hlt_tutorial/Muon2023C/93c501c3-7f03-4c42-806d-2e5c05d44364.root",
    "/eos/cms/store/group/dpg_trigger/comm_trigger/TriggerStudiesGroup/STEAM/hlt_tutorial/Muon2023C/93f39670-5c06-46cd-bedb-bf9217aec50b.root",
    "/eos/cms/store/group/dpg_trigger/comm_trigger/TriggerStudiesGroup/STEAM/hlt_tutorial/Muon2023C/96158cff-3839-41d6-962c-67d10d4df518.root",
    "/eos/cms/store/group/dpg_trigger/comm_trigger/TriggerStudiesGroup/STEAM/hlt_tutorial/Muon2023C/98bac889-850e-4b7d-a19f-22642b001903.root",
    "/eos/cms/store/group/dpg_trigger/comm_trigger/TriggerStudiesGroup/STEAM/hlt_tutorial/Muon2023C/9d91d42b-6e65-47a0-8ce8-304ab6bd2939.root",
    "/eos/cms/store/group/dpg_trigger/comm_trigger/TriggerStudiesGroup/STEAM/hlt_tutorial/Muon2023C/e26ad17a-14c8-409d-883b-7dbb16e65702.root",
    "/eos/cms/store/group/dpg_trigger/comm_trigger/TriggerStudiesGroup/STEAM/hlt_tutorial/Muon2023C/12793b56-216c-4119-a341-7477090e95fd.root",
    "/eos/cms/store/group/dpg_trigger/comm_trigger/TriggerStudiesGroup/STEAM/hlt_tutorial/Muon2023C/1891ca0e-ae5b-4272-9750-59fd1200fb3f.root",
    "/eos/cms/store/group/dpg_trigger/comm_trigger/TriggerStudiesGroup/STEAM/hlt_tutorial/Muon2023C/1eb1fb50-e4a4-4d2a-99ec-f33fa134cd28.root",
    "/eos/cms/store/group/dpg_trigger/comm_trigger/TriggerStudiesGroup/STEAM/hlt_tutorial/Muon2023C/22dc2eba-3c98-429a-8ae2-a22460d18e46.root",
    "/eos/cms/store/group/dpg_trigger/comm_trigger/TriggerStudiesGroup/STEAM/hlt_tutorial/Muon2023C/233e0089-6c14-48ab-b940-3caacab875c7.root",
    "/eos/cms/store/group/dpg_trigger/comm_trigger/TriggerStudiesGroup/STEAM/hlt_tutorial/Muon2023C/33f88913-81b6-4c08-8ffc-4b3c6d1c457d.root",
    "/eos/cms/store/group/dpg_trigger/comm_trigger/TriggerStudiesGroup/STEAM/hlt_tutorial/Muon2023C/5136f404-7feb-4381-b70e-c29c21e50063.root",
    "/eos/cms/store/group/dpg_trigger/comm_trigger/TriggerStudiesGroup/STEAM/hlt_tutorial/Muon2023C/60a36b0f-ba33-4946-b5c6-65eed9ff9d4a.root",
    "/eos/cms/store/group/dpg_trigger/comm_trigger/TriggerStudiesGroup/STEAM/hlt_tutorial/Muon2023C/67a4ae4c-d5bc-452d-baf0-66b6682f7d3c.root",
    "/eos/cms/store/group/dpg_trigger/comm_trigger/TriggerStudiesGroup/STEAM/hlt_tutorial/Muon2023C/7812f9a0-5122-476f-8cc6-c0f3e3de57d0.root",
    "/eos/cms/store/group/dpg_trigger/comm_trigger/TriggerStudiesGroup/STEAM/hlt_tutorial/Muon2023C/79d53f50-a484-44f0-8958-ddac1f07551e.root",
    "/eos/cms/store/group/dpg_trigger/comm_trigger/TriggerStudiesGroup/STEAM/hlt_tutorial/Muon2023C/7ed40a7c-737f-42bc-839f-cef8f5697d94.root",
    "/eos/cms/store/group/dpg_trigger/comm_trigger/TriggerStudiesGroup/STEAM/hlt_tutorial/Muon2023C/82e83a02-0a73-4632-b40b-317e4fa2cebe.root",
    "/eos/cms/store/group/dpg_trigger/comm_trigger/TriggerStudiesGroup/STEAM/hlt_tutorial/Muon2023C/8ecce0e1-5ca9-4c7a-a1e5-2552143b683b.root",
    "/eos/cms/store/group/dpg_trigger/comm_trigger/TriggerStudiesGroup/STEAM/hlt_tutorial/Muon2023C/962cb338-3e7c-4f5d-9f36-ce9f98cd9667.root",
    "/eos/cms/store/group/dpg_trigger/comm_trigger/TriggerStudiesGroup/STEAM/hlt_tutorial/Muon2023C/9a3acafe-fb88-4868-91a5-9c867e704ac6.root",
    "/eos/cms/store/group/dpg_trigger/comm_trigger/TriggerStudiesGroup/STEAM/hlt_tutorial/Muon2023C/a0c0ea1d-8d20-44d7-b63e-de92f732b595.root",
    "/eos/cms/store/group/dpg_trigger/comm_trigger/TriggerStudiesGroup/STEAM/hlt_tutorial/Muon2023C/a185faa2-23f8-4982-9812-94f394b6d84b.root",
    "/eos/cms/store/group/dpg_trigger/comm_trigger/TriggerStudiesGroup/STEAM/hlt_tutorial/Muon2023C/a47f9231-135a-4674-95a1-0b5647307b1a.root",
    "/eos/cms/store/group/dpg_trigger/comm_trigger/TriggerStudiesGroup/STEAM/hlt_tutorial/Muon2023C/a8812f4c-e695-4ccc-a8a3-ff00edd8610b.root",
    "/eos/cms/store/group/dpg_trigger/comm_trigger/TriggerStudiesGroup/STEAM/hlt_tutorial/Muon2023C/abea8734-580e-411a-ba3f-e0985bd73340.root",
    "/eos/cms/store/group/dpg_trigger/comm_trigger/TriggerStudiesGroup/STEAM/hlt_tutorial/Muon2023C/af3dd69f-7398-445b-a5cb-52a40082c97f.root",
    "/eos/cms/store/group/dpg_trigger/comm_trigger/TriggerStudiesGroup/STEAM/hlt_tutorial/Muon2023C/b4db366c-a3c6-48da-866d-0da53f0e8c53.root",
    "/eos/cms/store/group/dpg_trigger/comm_trigger/TriggerStudiesGroup/STEAM/hlt_tutorial/Muon2023C/c4d10eb4-8eed-404e-ab14-d8694fd5e785.root",
    "/eos/cms/store/group/dpg_trigger/comm_trigger/TriggerStudiesGroup/STEAM/hlt_tutorial/Muon2023C/c6e0037f-2d9b-49d3-a3c4-30df4264906c.root",
    "/eos/cms/store/group/dpg_trigger/comm_trigger/TriggerStudiesGroup/STEAM/hlt_tutorial/Muon2023C/c8862cec-5f41-4699-ba96-71e0b2736b30.root",
    "/eos/cms/store/group/dpg_trigger/comm_trigger/TriggerStudiesGroup/STEAM/hlt_tutorial/Muon2023C/c9a5b906-5615-4e1f-a447-41c2b8df0470.root",
    "/eos/cms/store/group/dpg_trigger/comm_trigger/TriggerStudiesGroup/STEAM/hlt_tutorial/Muon2023C/c9d8b7f5-a4d1-4c2b-8ca9-c78c26451eb8.root",
    "/eos/cms/store/group/dpg_trigger/comm_trigger/TriggerStudiesGroup/STEAM/hlt_tutorial/Muon2023C/dbc9f268-3288-4641-83ac-c5fddf0a28f3.root",
    "/eos/cms/store/group/dpg_trigger/comm_trigger/TriggerStudiesGroup/STEAM/hlt_tutorial/Muon2023C/e9fc4dda-f3dc-4d28-b93f-a6f88c06d93b.root",
    "/eos/cms/store/group/dpg_trigger/comm_trigger/TriggerStudiesGroup/STEAM/hlt_tutorial/Muon2023C/eb629065-21f7-420f-a6b5-c81621563d18.root",
    "/eos/cms/store/group/dpg_trigger/comm_trigger/TriggerStudiesGroup/STEAM/hlt_tutorial/Muon2023C/f1c02326-8cb3-441d-9f34-918fb1abc901.root",
    "/eos/cms/store/group/dpg_trigger/comm_trigger/TriggerStudiesGroup/STEAM/hlt_tutorial/Muon2023C/133d22fc-836d-4ac4-a28e-24cb4aa7461b.root",
    "/eos/cms/store/group/dpg_trigger/comm_trigger/TriggerStudiesGroup/STEAM/hlt_tutorial/Muon2023C/15b2e7c3-a9b4-4298-941e-c209a4a62d9e.root",
    "/eos/cms/store/group/dpg_trigger/comm_trigger/TriggerStudiesGroup/STEAM/hlt_tutorial/Muon2023C/1a0b8ccb-d1fb-4ba9-b4d6-ce38acffe54e.root",
    "/eos/cms/store/group/dpg_trigger/comm_trigger/TriggerStudiesGroup/STEAM/hlt_tutorial/Muon2023C/1cec4f58-70c9-4714-9a8e-160f25f0232a.root",
    "/eos/cms/store/group/dpg_trigger/comm_trigger/TriggerStudiesGroup/STEAM/hlt_tutorial/Muon2023C/21b53859-70b5-43f9-9ffd-f4e2ef3187a2.root",
    "/eos/cms/store/group/dpg_trigger/comm_trigger/TriggerStudiesGroup/STEAM/hlt_tutorial/Muon2023C/333457df-1d76-479e-9b97-0951f0808e77.root",
    "/eos/cms/store/group/dpg_trigger/comm_trigger/TriggerStudiesGroup/STEAM/hlt_tutorial/Muon2023C/35338b68-e715-4ae5-ad2d-8a8824119c32.root",
    "/eos/cms/store/group/dpg_trigger/comm_trigger/TriggerStudiesGroup/STEAM/hlt_tutorial/Muon2023C/3b502432-4343-491f-9da7-e0f4c28b9399.root",
    "/eos/cms/store/group/dpg_trigger/comm_trigger/TriggerStudiesGroup/STEAM/hlt_tutorial/Muon2023C/4930a204-66fe-4811-96ea-3d764852442f.root",
    "/eos/cms/store/group/dpg_trigger/comm_trigger/TriggerStudiesGroup/STEAM/hlt_tutorial/Muon2023C/4d2e67f5-1221-40c3-a068-fc2399821bac.root",
    "/eos/cms/store/group/dpg_trigger/comm_trigger/TriggerStudiesGroup/STEAM/hlt_tutorial/Muon2023C/4dd0acb6-7116-4cce-8223-dfeef4b0b14f.root",
    "/eos/cms/store/group/dpg_trigger/comm_trigger/TriggerStudiesGroup/STEAM/hlt_tutorial/Muon2023C/576080f2-39be-4f34-8d62-3d55ab9c5935.root",
    "/eos/cms/store/group/dpg_trigger/comm_trigger/TriggerStudiesGroup/STEAM/hlt_tutorial/Muon2023C/5e7b50aa-a169-4461-b219-3a03cb50aef5.root",
    "/eos/cms/store/group/dpg_trigger/comm_trigger/TriggerStudiesGroup/STEAM/hlt_tutorial/Muon2023C/6b40a9ae-5512-4205-b4c7-64ce5924c030.root",
    "/eos/cms/store/group/dpg_trigger/comm_trigger/TriggerStudiesGroup/STEAM/hlt_tutorial/Muon2023C/6b694fbf-3c6f-47f1-ad3c-24e999dc639d.root",
    "/eos/cms/store/group/dpg_trigger/comm_trigger/TriggerStudiesGroup/STEAM/hlt_tutorial/Muon2023C/75659266-dff1-4234-a641-8d1d3fa7e31a.root",
    "/eos/cms/store/group/dpg_trigger/comm_trigger/TriggerStudiesGroup/STEAM/hlt_tutorial/Muon2023C/7676b51f-aad0-4843-9195-47db42f551ac.root",
    "/eos/cms/store/group/dpg_trigger/comm_trigger/TriggerStudiesGroup/STEAM/hlt_tutorial/Muon2023C/77f3e032-9045-4be5-ace4-b95260ccc128.root",
    "/eos/cms/store/group/dpg_trigger/comm_trigger/TriggerStudiesGroup/STEAM/hlt_tutorial/Muon2023C/792946a6-dbd9-4512-a862-c81ad8167d8d.root",
    "/eos/cms/store/group/dpg_trigger/comm_trigger/TriggerStudiesGroup/STEAM/hlt_tutorial/Muon2023C/7b17b56d-17fb-445c-983b-d68268a0ee5e.root",
    "/eos/cms/store/group/dpg_trigger/comm_trigger/TriggerStudiesGroup/STEAM/hlt_tutorial/Muon2023C/7d74a13e-b680-44c9-bea9-7e91b96ea472.root",
]

reference_paths = ["Mu50", "IsoMu24"]

signal_paths    = ["AK8PFJet420_MassSD30",
                   "AK8PFJet425_SoftDropMass40",
                   "AK8PFJet450",
                   "AK8PFJet500",
                   "PFJet500",
                   "PFHT1050",
                   ]

p=PostProcessor(".",files,cut=preselection,branchsel=None,modules=[TrigHadAnalysis()],noOut=True,histFileName="histos_HadTrigNanoAOD.root",histDirName="hadTrigAnalyzerNanoAOD")
p.run()
