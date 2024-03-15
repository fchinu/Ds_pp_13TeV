from hipe4ml_converter.h4ml_converter import H4MLConverter
from hipe4ml.model_handler import ModelHandler
import argparse
import os

modelHandlers = [   '/home/fchinu/Run3/Ds_pp_13TeV/ML/Training/pt1_1.5/ModelHandler_pT_1_1.5.pickle',
                    '/home/fchinu/Run3/Ds_pp_13TeV/ML/Training/pt1.5_2/ModelHandler_pT_1.5_2.pickle',
                    '/home/fchinu/Run3/Ds_pp_13TeV/ML/Training/pt2_2.5/ModelHandler_pT_2_2.5.pickle',
                    '/home/fchinu/Run3/Ds_pp_13TeV/ML/Training/pt2.5_3/ModelHandler_pT_2.5_3.pickle',
                    '/home/fchinu/Run3/Ds_pp_13TeV/ML/Training/pt3_3.5/ModelHandler_pT_3_3.5.pickle',
                    '/home/fchinu/Run3/Ds_pp_13TeV/ML/Training/pt3.5_4/ModelHandler_pT_3.5_4.pickle',
                    '/home/fchinu/Run3/Ds_pp_13TeV/ML/Training/pt4_4.5/ModelHandler_pT_4_4.5.pickle',
                    '/home/fchinu/Run3/Ds_pp_13TeV/ML/Training/pt4.5_5/ModelHandler_pT_4.5_5.pickle',
                    '/home/fchinu/Run3/Ds_pp_13TeV/ML/Training/pt5_5.5/ModelHandler_pT_5_5.5.pickle',
                    '/home/fchinu/Run3/Ds_pp_13TeV/ML/Training/pt5.5_6/ModelHandler_pT_5.5_6.pickle',
                    '/home/fchinu/Run3/Ds_pp_13TeV/ML/Training/pt6_8/ModelHandler_pT_6_8.pickle',
                    '/home/fchinu/Run3/Ds_pp_13TeV/ML/Training/pt8_12/ModelHandler_pT_8_12.pickle',
                    '/home/fchinu/Run3/Ds_pp_13TeV/ML/Training/pt12_24/ModelHandler_pT_12_24.pickle'
                  ] 

# Get the model handlers
for modelPath in modelHandlers:
    model = ModelHandler()
    model.load_model_handler(modelPath)
    converter = H4MLConverter(model)
    converter.convert_model_onnx(1)
    converter.dump_model_onnx(modelPath.replace("pickle","onnx"))
    print(f"Model {modelPath} converted to ONNX")
