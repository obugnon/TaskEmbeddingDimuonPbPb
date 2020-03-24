/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef TaskDimuonPbPb_H
#define TaskDimuonPbPb_H

#include "AliAnalysisTaskSE.h"

class TaskDimuonPbPb : public AliAnalysisTaskSE  
{
    public:
                                TaskDimuonPbPb();
                                TaskDimuonPbPb(const char *name,int firstRun, int lastRun, UInt_t triggerClass);
        virtual                 ~TaskDimuonPbPb();
        
        virtual void            NotifyRun();
        Bool_t                  IsMuonFromJPsi(int muonIndex);
        virtual void            UserCreateOutputObjects();
        virtual void            UserExec(Option_t* option);
        virtual void            Terminate(Option_t* option);

    private:
        AliAODEvent*            fAODEvent;      //! input event
        AliVEvent*              fVEvent;        //! input event
        AliMuonTrackCuts*       fMuonTrackCuts; //! usual cuts on single muon tracks
        TClonesArray*           fArrayMCParticles;
        UInt_t                  fTriggerClass;  // trigger selection
        int                     fFirstRun, fLastRun;

        TH1I *fHistoTotalEventsPerRun;      // histogram to store number of events
        TH1I *fHistoPSEventsPerRun;         // for the physics selection
        TH1I *fHistoEventsBeforePSPerRun;   
      
        THnSparseD *fHistoDiMuonOS;     // histogram to store some properties of dimuons unlike signe      
        THnSparseD *fHistoSingleMuon;   // histogram to store some properties of single muons
        THnSparseD *fHistoJPsiGenerated;
        THnSparseD *fHistoSingleMuonGenerated;   

      
        TList *fListEventHistos;        // list to save the events histograms 
        TList *fListGeneratedHistos;
        TList *fListReconstructedHistos;

        TaskDimuonPbPb(const TaskDimuonPbPb&); // not implemented
        TaskDimuonPbPb& operator=(const TaskDimuonPbPb&); // not implemented

        ClassDef(TaskDimuonPbPb, 2);
};

#endif
