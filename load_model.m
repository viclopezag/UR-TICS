function model = load_model(metabolic_model)
% Load a genome scale model for the selected organism.
% Options:
% metabolic_model: 'iCG760', ' sMtb2.0
%                  'iSM810', 'sMtb',
%                  'GSMN_TB_1.1', 'iNJ661v_modified', 'iOSDD890'
%
% Author: Victor Lopez, Diciembre 10, 2018.


%% MODEL PATH

    VASHISHT_MODEL = 'models/iOSDD890.mat';
    GARAY_MODEL = 'models/iCG760_grRuleCorrected.mat'; 
    MA_MODEL = 'models/iSM810_grRuleCorrected.mat';
    LOFTHOUSE_MODEL = 'models/GSMN_TB1.1_grRuleCorrected.mat';
    RIENKSMA_MODEL = 'models/sMtb.mat';
    SMTB2_0_MODEL = 'models/sMtb2.0.mat';
    XAVIER_MODEL = 'models/iNJ661v_modified.mat';
    RIENKSMA_2018 = 'models/sMtb2018_EX.mat';
    KAVVAS_MODEL = 'models/iEK1011_biggdb.mat';
%% MEDIUM PATH

    MA_MEDIUM = 'medium/ExchangeRxns_MA2015.mat';
    LOFTHOUSE_MEDIUM = 'medium/ExchangeRxns_MA2015.mat';   

%% CHOICES    
    
switch metabolic_model
    
            case 'sMtb2.0'

                    load(SMTB2_0_MODEL);
                    model = sMtb2_0;
                    model.lb(model.lb < 0) = -1000;
                    model.ub(model.ub > 1) = 1000;
                    exchangeRxns  = model.rxns(findExcRxns(model)); 
                    model = changeRxnBounds(model,exchangeRxns,-1,'l'); 
                    model = changeRxnBounds(model,'Maintenance',1,'b'); % Maintenance reaction was fix to 1 mmol/gDW/h.
                    
          
					
            case 'iCG760'
                    
                    load(GARAY_MODEL);
                    model = iCG760;
                    model.lb(model.lb < 0) = -1000;
                    model.ub(model.ub > 1) = 1000;
                    exchangeRxns  = model.rxns(findExcRxns(model)); 
                    model = changeRxnBounds(model,exchangeRxns,-1,'l'); 
                    model = changeRxnBounds(model,'R129',1,'b'); % Mainteinence was fixt to 1 mmol/gDW/h
                    
    
            case 'iSM810'
                               
                    load (MA_MODEL);
                    model = iSM810;
                    model.lb(model.lb < 0) = -1000;
                    model.ub(model.ub > 1) = 1000; 
                    load (MA_MEDIUM);
                    model = changeRxnBounds(model,exchangeRxns,1,'u');
                    model = changeRxnBounds(model,'R129',1,'b'); % Maintainance was fix to 1 mmol/gDWH
                    
                    
                   
              case 'sMtb'
            
                    load(RIENKSMA_MODEL);
                    model = sMtb;
                    model.lb(model.lb < 0) = -1000;
                    model.ub(model.ub > 1) = 1000;
                    exchangeRxns  = model.rxns(findExcRxns(model)); 
                    model = changeRxnBounds(model,exchangeRxns,-1,'l'); 
                    model = changeRxnBounds(model,'Maintenance',1,'b'); % Maintenance reaction was fix to 1 mmol/gDW/h.
                    
                     
              case 'GSMN_TB_1.1'
                    
                    load(LOFTHOUSE_MODEL);
                    load(LOFTHOUSE_MEDIUM);
                    model = GSMN_TB1_1;
                    model.ub(model.ub > 1) = 1000;
                    model = changeRxnBounds(model,exchangeRxns,1,'u');
                    model = changeRxnBounds(model,'R129',1,'b'); % maintaince reaction was fix to 1 mmol/gDW/h

                    
              case 'iOSDD890'
                    
                    load(VASHISHT_MODEL);
                    model = iOSDD890;
                    model.lb(model.lb < 0) = -1000;
                    model.ub(model.ub > 1) = 1000;
                    exchangeRxns  = model.rxns(findExcRxns(model)); 
                    model = changeRxnBounds(model,exchangeRxns,-1,'l'); 
                    model = changeRxnBounds(model,'ATPS1',1,'b'); % maintaince reaction was fix to 1 mmol/gDW/h
					
                                 
               case 'iNJ661v_modified'
            
                    load(XAVIER_MODEL);
                    model = iNJ661v_modified;
                    model.lb(model.lb < 0) = -1000;
                    model.ub(model.ub > 1) = 1000;
                    exchangeRxns  = model.rxns(findExcRxns(model)); 
                    model = changeRxnBounds(model,exchangeRxns,-1,'l'); 
                    model = changeRxnBounds(model,'ATPS1',1,'b'); % maintaince reaction was fix to 1 mmol/gDW/h

               case 'iEK1011'
                    
                    load(KAVVAS_MODEL);
                    model = iEK1011;
                    model.lb(model.lb < 0) = -1000;
                    model.ub(model.ub > 1) = 1000;
                    exchangeRxns = model.rxns(findExcRxns(model));
                    model = changeRxnBounds(model,exchangeRxns,-1,'l');
                    model = changeRxnBounds(model,'ATPM',1,'b'); % maintaince reaction was fix to 1 mmol/gDW/h
              
                case 'sMtb2018'
                    
                    load(RIENKSMA_2018);
                    model = sMtb2018;
                    model.lb(model.lb < 0) = -1000;
                    model.ub(model.ub > 1) = 1000;
                    exchangeRxns  = model.rxns(findExcRxns(model)); 
                    model = changeRxnBounds(model,exchangeRxns,-1,'l'); 
                    model = changeRxnBounds(model,'R_Maintenance',1,'b'); % Maintenance reaction was fix to 1 mmol/gDW/h.
                    
                  
end
              
              
end