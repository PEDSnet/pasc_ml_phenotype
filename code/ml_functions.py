from sklearn.metrics import precision_recall_curve, roc_auc_score, f1_score, make_scorer, recall_score, accuracy_score, precision_score, confusion_matrix, average_precision_score
import pandas as pd
import numpy as np

def model_eval(model, X_train, y_train, param_grid={}, 
               scoring={'accuracy': make_scorer(accuracy_score), 'precision': make_scorer(precision_score),
         'recall': make_scorer(recall_score),'roc_auc': make_scorer(roc_auc_score, needs_proba=True),
         'average_precision': make_scorer(average_precision_score), 'f1': make_scorer(f1_score)},
               refit="roc_auc"):
    if param_grid=={}:
        return pd.DataFrame(cross_validate(model, X_train, y_train, scoring=scoring))#.mean()
    else:
        fit_grid=GridSearchCV(model, param_grid, scoring=scoring, refit=refit, verbose=1).fit(X_train, y_train)
        return pd.DataFrame(fit_grid.cv_results_)
def model_score(model, X_train, y_train, X_test, y_test):
    fit=model.fit(X_train, y_train)
    y_pred=model.predict(X_test)
    y_prob=model.predict_proba(X_test)[:,1]
    acc=accuracy_score(y_test, y_pred)
    rec=recall_score(y_test, y_pred)
    prec=precision_score(y_test, y_pred)
    auroc=roc_auc_score(y_test, y_prob)
    ap=average_precision_score(y_test, y_prob)
    f1=f1_score(y_test, y_pred)
    score={'accuracy':acc, 'recall':rec, 'precision':prec, 'auroc':auroc, 'average precision':ap, 'f1':f1}
    return score

def score_nn(trained_model, X_input, y_test, thresh):
    y_prob=trained_model.predict(X_input)
    y_pred=y_prob.copy()
    y_pred[y_prob<thresh]=0
    y_pred[y_prob>=thresh]=1
    acc=accuracy_score(y_test, y_pred)
    rec=recall_score(y_test, y_pred)
    prec=precision_score(y_test, y_pred)
    auroc=roc_auc_score(y_test, y_prob)
    ap=average_precision_score(y_test, y_prob)
    f1=f1_score(y_test, y_pred)
    score={'accuracy':acc, 'recall':rec, 'precision':prec, 'auroc':auroc, 'average precision':ap, 'f1':f1}
    return score

def score_xg_multioutput(trained_model, X_input, y_test, thresh):
    y_test_pasc=y_test[:,0]
    y_test_non_misc=y_test[:,1]
    y_prob=trained_model.predict_proba(X_input)
    y_prob_pasc=y_prob[0][:,1]
    y_prob_non_misc=y_prob[1][:,1]
    y_pred_pasc=y_prob_pasc.copy()
    y_pred_non_misc=y_prob_non_misc.copy()
    y_pred_pasc[y_prob_pasc<thresh]=0
    y_pred_pasc[y_prob_pasc>=thresh]=1
    y_pred_non_misc[y_prob_non_misc<thresh]=0
    y_pred_non_misc[y_prob_non_misc>=thresh]=1
    acc_pasc=accuracy_score(y_test_pasc, y_pred_pasc)
    rec_pasc=recall_score(y_test_pasc, y_pred_pasc)
    prec_pasc=precision_score(y_test_pasc, y_pred_pasc)
    auroc_pasc=roc_auc_score(y_test_pasc, y_prob_pasc)
    ap_pasc=average_precision_score(y_test_pasc, y_prob_pasc)
    f1_pasc=f1_score(y_test_pasc, y_pred_pasc)
    
    acc_non_misc=accuracy_score(y_test_non_misc, y_pred_non_misc)
    rec_non_misc=recall_score(y_test_non_misc, y_pred_non_misc)
    prec_non_misc=precision_score(y_test_non_misc, y_pred_non_misc)
    auroc_non_misc=roc_auc_score(y_test_non_misc, y_prob_non_misc)
    ap_non_misc=average_precision_score(y_test_non_misc, y_prob_non_misc)
    f1_non_misc=f1_score(y_test_non_misc, y_pred_non_misc)

    score={'accuracy_pasc':acc_pasc, 'recall_pasc':rec_pasc, 'precision_pasc':prec_pasc,
                'auroc_pasc':auroc_pasc, 'average precision_pasc':ap_pasc, 'f1_pasc':f1_pasc,
               'accuracy_non_misc':acc_non_misc, 'recall_non_misc':rec_non_misc, 'precision_non_misc':prec_non_misc,
                'auroc_non_misc':auroc_non_misc, 'average precision_non_misc':ap_non_misc, 'f1_non_misc':f1_non_misc}
    return score

def get_scores(y_prob_pasc_any, y_prob_misc, y_prob_non_misc, y_test_pasc_any, y_test_misc, y_test_non_misc, thresh):
    y_pred_non_misc=(y_prob_non_misc>thresh)*1
    y_pred_misc=(y_prob_misc>thresh)*1
    y_pred_pasc_any=(y_prob_pasc_any>thresh)*1
    
    acc_pasc_any=accuracy_score(y_test_pasc_any, y_pred_pasc_any)
    rec_pasc_any=recall_score(y_test_pasc_any, y_pred_pasc_any)
    spec_pasc_any=recall_score(y_test_pasc_any, y_pred_pasc_any, pos_label=0)
    prec_pasc_any=precision_score(y_test_pasc_any, y_pred_pasc_any)
    auroc_pasc_any=roc_auc_score(y_test_pasc_any, y_prob_pasc_any)
    ap_pasc_any=average_precision_score(y_test_pasc_any, y_prob_pasc_any)
    f1_pasc_any=f1_score(y_test_pasc_any, y_pred_pasc_any)
    youden_j_pasc_any=rec_pasc_any+spec_pasc_any-1
    
    acc_misc=accuracy_score(y_test_misc, y_pred_misc)
    rec_misc=recall_score(y_test_misc, y_pred_misc)
    spec_misc=recall_score(y_test_misc, y_pred_misc, pos_label=0)    
    prec_misc=precision_score(y_test_misc, y_pred_misc)
    auroc_misc=roc_auc_score(y_test_misc, y_prob_misc)
    ap_misc=average_precision_score(y_test_misc, y_prob_misc)
    f1_misc=f1_score(y_test_misc, y_pred_misc)
    youden_j_misc=rec_misc+spec_misc-1

    
    acc_non_misc=accuracy_score(y_test_non_misc, y_pred_non_misc)
    rec_non_misc=recall_score(y_test_non_misc, y_pred_non_misc)
    spec_non_misc=recall_score(y_test_non_misc, y_pred_non_misc, pos_label=0)
    prec_non_misc=precision_score(y_test_non_misc, y_pred_non_misc)
    auroc_non_misc=roc_auc_score(y_test_non_misc, y_prob_non_misc)
    ap_non_misc=average_precision_score(y_test_non_misc, y_prob_non_misc)
    f1_non_misc=f1_score(y_test_non_misc, y_pred_non_misc)
    youden_j_non_misc=rec_non_misc+spec_non_misc-1


    score={'accuracy_pasc_any':acc_pasc_any, 'recall_pasc_any':rec_pasc_any,
           'specificity_pasc_any':spec_pasc_any, 'precision_pasc_any':prec_pasc_any,
           'auroc_pasc_any':auroc_pasc_any, 'average precision_pasc_any':ap_pasc_any, 'f1_pasc_any':f1_pasc_any,
           'youden_j_pasc_any':youden_j_pasc_any,
           'accuracy_non_misc':acc_non_misc, 'recall_non_misc':rec_non_misc, 
           'specificity_non_misc':spec_non_misc, 'precision_non_misc':prec_non_misc,
           'auroc_non_misc':auroc_non_misc, 'average precision_non_misc':ap_non_misc, 'f1_non_misc':f1_non_misc,
           'youden_j_non_misc':youden_j_non_misc,
           'accuracy_misc':acc_misc, 'recall_misc':rec_misc, 
           'spec_misc':spec_misc, 'precision_misc':prec_misc,
           'auroc_misc':auroc_misc, 'average precision_misc':ap_misc, 'f1_misc':f1_misc,
          'youden_j_misc':youden_j_misc}
    return score

def youden_j_pasc_any(results_list, p):
    scores=get_scores(results_list[0][2]['y_prob_pasc_any'], results_list[0][2]['y_prob_misc'], 
           results_list[0][2]['y_prob_non_misc'], results_list[0][1]['y_test_pasc_any'], 
           results_list[0][1]['y_test_misc'], results_list[0][1]['y_test_non_misc'], thresh=p)
    return(scores['youden_j_pasc_any'])


def score_xg_multiclass(trained_model, X_input, y_test_multiclass, thresh):
    y_test_non_misc=(y_test_multiclass==2)*1
    y_test_misc=(y_test_multiclass==1)*1
    y_test_pasc_any=((y_test_multiclass==1) | (y_test_multiclass==2))*1
    y_prob=trained_model.predict_proba(X_input)
    y_prob_non_misc=y_prob[:,2]
    y_prob_misc=y_prob[:,1]
    y_prob_pasc_any=y_prob[:,1]+y_prob[:,2]
    y_pred_non_misc=(y_prob_non_misc>thresh)*1
    y_pred_misc=(y_prob_misc>thresh)*1
    y_pred_pasc_any=(y_prob_pasc_any>thresh)*1
    
    acc_pasc_any=accuracy_score(y_test_pasc_any, y_pred_pasc_any)
    rec_pasc_any=recall_score(y_test_pasc_any, y_pred_pasc_any)
    prec_pasc_any=precision_score(y_test_pasc_any, y_pred_pasc_any)
    auroc_pasc_any=roc_auc_score(y_test_pasc_any, y_prob_pasc_any)
    ap_pasc_any=average_precision_score(y_test_pasc_any, y_prob_pasc_any)
    f1_pasc_any=f1_score(y_test_pasc_any, y_pred_pasc_any)
    
    acc_misc=accuracy_score(y_test_misc, y_pred_misc)
    rec_misc=recall_score(y_test_misc, y_pred_misc)
    prec_misc=precision_score(y_test_misc, y_pred_misc)
    auroc_misc=roc_auc_score(y_test_misc, y_prob_misc)
    ap_misc=average_precision_score(y_test_misc, y_prob_misc)
    f1_misc=f1_score(y_test_misc, y_pred_misc)
    
    acc_non_misc=accuracy_score(y_test_non_misc, y_pred_non_misc)
    rec_non_misc=recall_score(y_test_non_misc, y_pred_non_misc)
    prec_non_misc=precision_score(y_test_non_misc, y_pred_non_misc)
    auroc_non_misc=roc_auc_score(y_test_non_misc, y_prob_non_misc)
    ap_non_misc=average_precision_score(y_test_non_misc, y_prob_non_misc)
    f1_non_misc=f1_score(y_test_non_misc, y_pred_non_misc)

    score={'accuracy_pasc_any':acc_pasc_any, 'recall_pasc_any':rec_pasc_any, 'precision_pasc_any':prec_pasc_any,
           'auroc_pasc_any':auroc_pasc_any, 'average precision_pasc_any':ap_pasc_any, 'f1_pasc_any':f1_pasc_any,
           'accuracy_non_misc':acc_non_misc, 'recall_non_misc':rec_non_misc, 'precision_non_misc':prec_non_misc,
           'auroc_non_misc':auroc_non_misc, 'average precision_non_misc':ap_non_misc, 'f1_non_misc':f1_non_misc,
           'accuracy_misc':acc_misc, 'recall_misc':rec_misc, 'precision_misc':prec_misc,
           'auroc_misc':auroc_misc, 'average precision_misc':ap_misc, 'f1_misc':f1_misc}
    y_test={'y_test_pasc_any':y_test_pasc_any, 'y_test_misc':y_test_misc, 'y_test_non_misc':y_test_non_misc}
    y_prob={'y_prob_pasc_any':y_prob_pasc_any, 'y_prob_misc':y_prob_misc, 'y_prob_non_misc':y_prob_non_misc}
    results=[score, y_test, y_prob]
    return results


def make_feature_df(demo, cond_features, drug_features, util_features, lab_features, proc_features, lab_cts):
    site=pd.get_dummies(demo.site, prefix='site')
    sex=pd.get_dummies(demo.sex_cat, prefix='sex')
    eth=pd.get_dummies(demo.eth_cat, prefix='eth')
    period=pd.get_dummies(demo.cohort_entry_period, prefix='cohort_entry')
    X_demo=pd.concat([demo.pasc_flag, demo.pasc_not_misc_flag, site, sex, eth, period, demo.entry_age], axis=1)
    X_temp=pd.merge(pd.merge(pd.merge(pd.merge(pd.merge(pd.merge(X_demo, cond_features, left_index=True, right_index=True, how='inner'),drug_features, left_index=True, right_index=True, how='inner'), util_features, left_index=True, right_index=True, how='inner'), lab_features, left_index=True, right_index=True, how='inner'), proc_features, left_index=True, right_index=True, how='inner'), lab_cts, left_index=True, right_index=True, how='inner')
    y_pasc=X_temp.pasc_flag.values
    y_pasc_not_misc=X_temp.pasc_not_misc_flag.values
    y=np.column_stack((y_pasc, y_pasc_not_misc))
    #y=y_pasc
    X=X_temp.drop(columns=['pasc_flag', 'pasc_not_misc_flag'])
    return X, y
def drop_diff(X_train, X_test):
    X_test=X_test.drop(columns=X_test.columns.difference(X_train.columns))
    X_train=X_train.drop(columns=X_train.columns.difference(X_test.columns))
    return X_train, X_test

def make_path(domain, train_or_test, i, path_start="../specs/cohort_features/"):
    if (domain=="demo"):
        if (i=="all"):
            return("".join([path_start, "cohort_untested_demo_", str(i), ".csv"]))
        else:
            return("".join([path_start, "cohort_untested_demo_", train_or_test, "_", str(i), ".csv"]))
    elif (domain=="lab_cts"):
        if (i=="all"):
            return("".join([path_start, "cohort_untested_lab_cts_", str(i), ".csv"]))
        else:
            return("".join([path_start, "cohort_lab_cts_", train_or_test, "_", str(i), ".csv"]))
    elif (domain=="util" or domain=="drug" or domain=="cond" or domain=="proc" or domain=="lab"):
        if (domain=="util" and i=="all"):
            return("".join([path_start, "cohort_untested_", domain,  "_",str(i), ".csv"]))
        elif (i=="all"):
            return("".join([path_start, "cohort_", domain, "_features_",  str(i), ".csv"]))
        else:
            return("".join([path_start, "cohort_", domain, "_features_", train_or_test, "_", str(i), ".csv"]))
    else:
        return("")
