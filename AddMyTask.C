#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TString.h"
#include "TObjArray.h"

#include "AliLog.h"
#include "AliVEventHandler.h"

#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"

#include "AliMuonTrackCuts.h"
#include "TaskDimuonPbPb.h"
#endif

TString GetTriggerClassName(UInt_t triggerClass = AliVEvent::kINT7inMUON);

///Add the names of the output containers. The method is not fully automatic with the cxx file, so you have to keep the same order here and there (when defining the output containers in TaskDimuonPbPb.cxx)
TString outputListNames[] = {"EventHistos", "SingleMuonHistos", "DiMuonHistos"};
int numberOfOutputLists = sizeof( outputListNames ) / sizeof( outputListNames[0] );

///Add the names of the trigger classes. Refer to AliVEvent.h to see the description for each trigger class. When adding multiple entries in the array, the analysis manager will run the class multiple times, each time using the correponding trigger class. The output containers will be labelled by the name of trigger
UInt_t triggerClasses[] = {AliVEvent::kAny};
int numberOfTriggerClasses = sizeof( triggerClasses ) / sizeof( triggerClasses[0] );
const int numberOfTriggerClasses_CONST = numberOfTriggerClasses;

TaskDimuonPbPb* AddMyTask(Bool_t isESD=kFALSE, TString name = "name", int firstRun = 244918, int lastRun = 297595)
{
    // get the manager via the static access member. since it's static, you don't need
    // to create an instance of the class here to call the function
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        Error("AddTaskDimuonPbPb","AliAnalysisManager not set!");
        return 0x0;
    }
    // get the input event handler, again via a static method. 
    // this handler is part of the managing system and feeds events to your task
    TString type = mgr->GetInputEventHandler()->GetDataType();
    if (!type.Contains("AOD") && !isESD) 
    {
        Error("AddTaskDimuonPbPb", "AOD input handler needed!");
        return 0x0;
    }
    if (!type.Contains("ESD") && isESD) 
    {
        Error("AddTaskDimuonPbPb", "ESD input handler needed!");
        return 0x0;
    }

    // by default, a file is open for writing. here, we get the filename
    TString outputFileName = AliAnalysisManager::GetCommonFileName();

    TaskDimuonPbPb *task[numberOfTriggerClasses_CONST];
    for(int iTask =0;iTask<numberOfTriggerClasses;iTask++)
    {
        TString triggerClassName = GetTriggerClassName(triggerClasses[iTask]);
        //A task for each trigger class
        task[iTask] = new TaskDimuonPbPb(Form("TaskDimuonPbPb%s",triggerClassName.Data()) , firstRun , lastRun , triggerClasses[iTask]);
        if(!task[iTask]) return 0x0;
        //task[iTask]->SelectCollisionCandidates(triggerClasses[iTask]);
        // add your task to the manager
        mgr->AddTask(task[iTask]);
        // your task needs input: here we connect the manager to your task
        mgr->ConnectInput(task[iTask], 0, mgr->GetCommonInputContainer());

        for(int iList =0;iList<numberOfOutputLists;iList++)
        {
            AliAnalysisDataContainer *outputList = mgr->CreateContainer(Form("%s_%s",outputListNames[iList].Data(),triggerClassName.Data()), TList::Class(), AliAnalysisManager::kOutputContainer,outputFileName);
            // same for the output
            mgr->ConnectOutput(task[iTask], iList+1, outputList);
        }
    }
    return task[0];
}

TString GetTriggerClassName(UInt_t triggerClass)
{

    //Dummy function to associate a name for each trigger event class. The name is to be used for naming the output container
    UInt_t availableTriggerClasses[] = {AliVEvent::kINT7inMUON, AliVEvent::kMuonSingleLowPt7, AliVEvent::kMuonSingleHighPt7, AliVEvent::kMuonLikeLowPt7, AliVEvent::kMuonUnlikeLowPt7, AliVEvent::kAny}; 

    int numberOfAvailableTriggerClasses = sizeof( availableTriggerClasses ) / sizeof( availableTriggerClasses[0] );

    TString availableTriggerClassesNames[] = {"CINT7", "CMSL7", "CMSH7", "CMLL7", "CMUL7", "CAny"};


    for(int iTrigger=0;iTrigger<numberOfAvailableTriggerClasses;iTrigger++){
      if(triggerClass == availableTriggerClasses[iTrigger]) return availableTriggerClassesNames[iTrigger];
    }

    return "";
}