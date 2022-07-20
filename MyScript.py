import ROOT
import helper
from helper import *
import csv
import json
import os
import pickle
import operator
from helper import deltaR , matching
from configure import *
from ROOT import TFile, TTree, gRandom
from array import array
from optparse import OptionParser

parser = OptionParser()
parser.add_option("--sname", dest="sname", default="G1Jet_Pt", action="store", help="can be QCD , GJets_Pt ... ")
parser.add_option("--stype", dest="stype", default="signal", action="store", help="can be data or signal or bkg")
parser.add_option("--ndiv", dest="ndiv", default="10", action="store", help="number of divitions for one root file")
parser.add_option("--divIndex", dest="divIndex", default="0", action="store", help="index of divitions for one root file")
(options, args) = parser.parse_args()
exec("ndiv="+options.ndiv)
exec("divIndex="+options.divIndex)
stype = options.stype
sname = options.sname

pfile = "samples_ana.pkl"
sample_dic = pickle.load(open(pfile,'rb'))
mychain_dict  =  getChain(year=2016, stype=stype, sname=sname, pfile=pfile, datatype='all', test=False)
ch = mychain_dict[0]

if stype == "signal" and sname == "G1Jet_Pt" :
	ch.Draw(">>eList", "ngoodPhoton==1&&(goodPhoton_pt>=225)&&(ngoodbJet==0)&&(abs(goodGenPhoton_pt-goodPhoton_pt)/goodPhoton_pt<0.1)") # Apply ROOT cuts here
 	elist = ROOT.gDirectory.Get("eList")
 	number_events = elist.GetN()
	#number_events = 1000
if stype == "bkg" and sname == "QCD_HT" :
        ch.Draw(">>eList", "ngoodPhoton==1&&(goodPhoton_pt>=225)&&(ngoodbJet==0)&&ngoodGenPhoton==0") # Apply ROOT cuts here
        elist = ROOT.gDirectory.Get("eList")
        number_events = elist.GetN()

nEventsPerChunk = number_events/float(ndiv)
ini_event = divIndex*int(nEventsPerChunk)
fin_event = min((divIndex+1)*int(nEventsPerChunk),number_events)

#Number of Particles in n events
Nph=set()  # N Photons per Event in goodPhotons
Nbj=set()
Nj=set()
NGbj=set()

#List of Variables for The dataset

NPhotons=[] # NPhotons before Matching
nPhotons=[] # nPhotons after Matching
nGPhotons=[] # nGenPhotons after Matching 
GPhotons=[] # List of GenPhotns variables
RGPhotons=[] # list of 'Raw' GenPhotons before Matching
RRPhotons=[] # List of Photons before matching
Sel_photons=[] # List of Photons after Matching

nJets=[] # nJets after selection
Sel_jets=[] # Lead Jets after selection
Sublead_jets=[]	# SubLead Jets after Selection

nbJets=[] # 1/2 nbJets before Matching 
Sel_bjets=[] # Lead bJets after Matching 
SubLead_bJets=[] # SubLead bJets after Matching
RLead_bJets=[] # Lead bJets before Matching 
RSubLead_bJets=[] # SubLead bJets before Matching

Gbjets=[] # Lead GenbJets after Matching
SubLead_GbJets=[] # SubLead GenbJets after Matching
RGbjets=[] # Lead GenbJets before Matching
RSubLead_GbJets=[] # SubLead GenbJets before Matching 

Met_pt=[] # Met_pt after Matching
Met_phi=[]
DRm=[]
Labels=[]
#Match_ph=[]
Match_1b=[] # 1b Jets after Matching
Match_2b=[] # 2b Jets after Matching
Number_b=0. # Number 1/2b before Matching
M_1b=0. # N 1b after
M_2b=0. # N 2b after
npho=0. # N Pho After 
#M_ph=0.
Npho=0. # N Pho  before
#Npho1=0.
#Npho2=0.
count_events_num=0
count_events_den=0
#--------- First loop : Looping on Events ---------

for jentry in range(ini_event,fin_event):
	#ch.GetEntry(jentry)
	ch.GetEntry(elist.GetEntry(jentry))

	good_event_Photon=False
	good_event_b=False
	good_event_2b=False

   	nGenPhoton=ch.GetLeaf('ngoodGenPhoton').GetValue() # nGen Particles in one event#float
   	nGenJet=ch.GetLeaf('nGenJet').GetValue() # nGen Particles in one event
   	nPhoton = ch.GetLeaf('ngoodPhoton').GetValue() # nReco Photon inside one event
   	nJet = ch.GetLeaf('ngoodJet').GetValue() # nReco Jets inside one event
   	nbJet = ch.GetLeaf('ngoodbJet').GetValue() # nReco Jets inside one event
   	met_pt=ch.GetLeaf('MET_pt').GetValue()
   	met_phi=ch.GetLeaf('MET_phi').GetValue()

	Nph.add(nPhoton)
   	Nbj.add(nbJet)
   	Nj.add(nJet)
	
	

	#if nPhoton != 1 or nGenPhoton  != 1 or nbJet != 0 : continue     # add control plotter cuts

	#count_events_den += 1 # All 1 Photon 1 /2 bJet Events

	#lists of variables per event
	Photons=[]          # kin of Reco Photons
   	Jets=[]             # kin of Reco Jets
   	GenJets=[]              # kin of G jets
	GenbJets=[]
   	GenPhotons=[]           # kin of Gen Particles
   	bJets=[]

	matched_photons=[]
   	sel_photons=[]
	matched_jets_1b=[]
	matched_jets_lb=[]
	matched_jets_slb=[]
	sel_jets_1b=[]
	sel_jets_lb=[]
	sel_jets_slb=[]
	lead_b=[]
	sublead_b=[]
	Glead_b=[]
	Lead_bJets=[]
	subGlead_b=[]
	dRms=[]
	labels=[]
	
	count_events_den += 1
       	nbJets.append({'nbJets':nbJet})
      	NPhotons.append(nPhoton)
	
#---------- second loop : nVariables inside Particles-------------------
	temp_index = 0
	
	for ID in range(int(nGenPhoton)): # Gen Photons   #List of Dict
			GenPhotons.append({'index':temp_index,'orig_index':ID,'pt':ch.GetLeaf('goodGenPhoton_pt').GetValue(ID),'eta':ch.GetLeaf('goodGenPhoton_eta').GetValue(ID),'phi':ch.GetLeaf('goodGenPhoton_phi').GetValue(ID)})
			temp_index += 1

	for I in range(int(nPhoton)): # Reco PHOTONS
			#Looping on ngoodPhotons Make sure that other Variables corresponds to a goodPhoton
			Photons.append({'index':I,'pt':ch.GetLeaf('goodPhoton_pt').GetValue(I),'eta':ch.GetLeaf('goodPhoton_eta').GetValue(I),'phi':ch.GetLeaf('goodPhoton_phi').GetValue(I),'hoe':ch.GetLeaf('Photon_hoe').GetValue(I),'sieie':ch.GetLeaf('Photon_sieie').GetValue(I),'r9':ch.GetLeaf('Photon_r9').GetValue(I),'Iso_all':ch.GetLeaf('Photon_pfRelIso03_all').GetValue(I),'Iso_chg':ch.GetLeaf('Photon_pfRelIso03_chg').GetValue(I)})

	for j in range(int(nJet)): # Reco JETS
        	Jets.append({'index':j,'pt':ch.GetLeaf('goodJet_pt').GetValue(j),'phi':ch.GetLeaf('goodJet_phi').GetValue(j),'eta':ch.GetLeaf('goodJet_eta').GetValue(j)})

	for j in range(int(nbJet)): # Reco  bJets
		bJets.append({'index':j,'pt':ch.GetLeaf('goodbJet_pt').GetValue(j),'phi':ch.GetLeaf('goodbJet_phi').GetValue(j),'eta':ch.GetLeaf('goodbJet_eta').GetValue(j)})


	print(' ')
	print('event#',jentry)
   	print('GenPhotons',len(GenPhotons),'Photons',len(Photons),'GenJets' , len(GenJets) ,'Jets',len(Jets),'Gen bJets',len(GenbJets),'bJets' ,len(bJets))
	print(' ')

# ------------PHOTON MATCHING INSIDE CONE dR----------------------------------------------
	
	RRPhotons += Photons
		
    	if  len(Photons)  :
            	sel_photons = Photons
            	good_event_Photon=True
		
	if good_event_Photon==True:
		labels.append({'label':'s'})
	else:
		labels.append({'label':'b'})
	print(labels)
	

#------------- Getting Leading & SubLeading Jet---------------------

    	lead_jets=sorted(Jets,key=lambda x:x['pt'],reverse=True)[:1]
	sublead_jets=sorted(Jets,key=lambda x:x['pt'],reverse=True)[1:2]
	
#-----------Extract Selected Variables----------------------------
	Match_1b.append(0)
        Match_2b.append(len(lead_b)+len(sublead_b))
        #Match_ph.append(len(sel_photons))	
        if len(matched_photons)==0:
                matched_photons.append({'pt':-999,'eta':-999,'phi':-999})
	GPhotons+=matched_photons
	if len(GenPhotons)==0:
              	RGPhotons.append({'pt':-999,'eta':-999,'phi':-999})
	RGPhotons+=GenPhotons
	if len(Photons)==0:
                RRPhotons.append({'pt':-999,'eta':-999,'phi':-999})
	#RRPhotons+=Photons		
        if len(GenPhotons)==0:
                nGPhotons.append({'nGPhotons':0})
	nGPhotons.append({'nGPhotons':len(GenPhotons)})
	if len(Photons)==0:
                nPhotons.append({'nPhotons':0})
        nPhotons.append({'nPhotons':len(Photons)})		
        if len(sel_photons)==0:
                sel_photons.append({'pt':-999,'eta':-999,'phi':-999,'hoe':-999,'sieie':-999,'r9':-999,'Iso_all':-999,'Iso_chg':-999})
	Sel_photons+=sel_photons		
        if len(Jets)==0:
                nJets.append({'nJets':0})
	nJets.append({'nJets':len(Jets)})				        
        if len(lead_jets)==0:
                lead_jets.append({'pt':-999,'eta':-999,'phi':-999})	
	Sel_jets+=lead_jets			        
        if len(sublead_jets)==0:
                sublead_jets.append({'pt':-9999,'eta':-9999,'phi':-9999})
	Sublead_jets+=sublead_jets 
        if met_pt==0:
                met_pt=-999						
        Met_pt.append({'pt':met_pt})
        if met_phi==0:
                met_phi=-999
        Met_phi.append({'phi':met_phi})
        if len(labels)==0:
                labels.append({'label':-999})
        Labels+=labels
	if good_event_Photon : # num cut satisfied
		count_events_num += 1
		
		

for n in nPhotons:
                npho += n['nPhotons']

#---------Create Dataset's list---------------------------------
written=[]
biggest_len = len(GPhotons)
for i in range(biggest_len):
	temp_dict = {} #empty dict
	temp_dict['RGPhoton_pt'] = RGPhotons[i]['pt']
        temp_dict['RGPhoton_eta'] = RGPhotons[i]['eta']
        temp_dict['RGPhoton_phi'] = RGPhotons[i]['phi']
        temp_dict['Rphoton_pt'] = RRPhotons[i]['pt']
        temp_dict['Rphoton_eta'] = RRPhotons[i]['eta']
        temp_dict['Rphoton_phi'] = RRPhotons[i]['phi']
        temp_dict['GPhoton_pt'] = GPhotons[i]['pt']
        temp_dict['GPhoton_eta'] = GPhotons[i]['eta']
        temp_dict['GPhoton_phi'] = GPhotons[i]['phi']
        temp_dict['photon_pt'] = Sel_photons[i]['pt']
        temp_dict['photon_eta'] = Sel_photons[i]['eta']
        temp_dict['photon_phi'] = Sel_photons[i]['phi']
	temp_dict['photon_hoe'] = Sel_photons[i]['hoe']
        temp_dict['photon_sieie'] = Sel_photons[i]['sieie']
        temp_dict['photon_r9'] = Sel_photons[i]['r9']
	temp_dict['photon_Iso_all'] = Sel_photons[i]['Iso_all']
        temp_dict['photon_Iso_chg'] = Sel_photons[i]['Iso_chg']
        temp_dict['jet_pt'] = Sel_jets[i]['pt']
        temp_dict['jet_eta'] = Sel_jets[i]['eta']
        temp_dict['jet_phi'] = Sel_jets[i]['phi']
        temp_dict['SLjet_pt'] = Sublead_jets[i]['pt']
        temp_dict['SLjet_eta'] = Sublead_jets[i]['eta']
        temp_dict['SLjet_phi'] = Sublead_jets[i]['phi']

        temp_dict['met_pt'] = Met_pt[i]['pt']
        temp_dict['met_phi'] = Met_phi[i]['phi']
        temp_dict['nGPhotons'] = nGPhotons[i]['nGPhotons']
        temp_dict['Labels'] = Labels[i]['label']	
	written.append(temp_dict)


#------- Write to CSV -------------------

fields = [ 'RMatched_GPhoton_pt','RMatched_GPhoton_eta','RMatched_GPhoton_phi','RPhoton_pt','RPhoton_eta','RPhoton_phi','Matched_GPhoton_pt','Matched_GPhoton_eta','Matched_GPhoton_phi','Photon_pt','Photon_eta','Photon_phi','Photon_hoe','Photon_sieie','Photon_r9','Iso_all','Iso_chg' ,'Lead_Jet_pt','Lead_Jet_eta','Lead_Jet_phi','SubLead_Jet_pt','SubLead_Jet_eta','SubLead_Jet_phi','MET_pt','MET_phi','Labels']


rows = written
with open('/afs/cern.ch/user/m/mbarakat/gamma_b/paraCSVQCD/Photon0b_Train'+str(divIndex)+'_'+str(stype)+'.csv', 'w') as f:

    		write = csv.writer(f)
    		write.writerow(fields)
    		for row in rows:
			r = [row['RGPhoton_pt'],row['RGPhoton_eta'],row['RGPhoton_phi'],row['Rphoton_pt'],row['Rphoton_eta'],row['Rphoton_phi'],row['GPhoton_pt'],row['GPhoton_eta'],row['GPhoton_phi'],row['photon_pt'],row['photon_eta'],row['photon_phi'],row['photon_hoe'],row['photon_sieie'],row['photon_r9'],row['photon_Iso_all'],row['photon_Iso_chg'], row['jet_pt'],row['jet_eta'],row['jet_phi'],row['SLjet_pt'],row['SLjet_eta'],row['SLjet_phi'],row['met_pt'],row['met_phi'],row['Labels']]
                        
        		write.writerow(r)




#------- Check for n objects in goodParticles--------
print('Nph',Nph)   # Check for extra photons in ngoodPhotons
print('nbJet', Nbj)
print('nGbj', NGbj)
print('count_events_num',count_events_num)
print('count_events_den',count_events_den)
Matching_Events_eff=(count_events_num)/(count_events_den*1.0)
print('Matching_Events_eff',Matching_Events_eff)
