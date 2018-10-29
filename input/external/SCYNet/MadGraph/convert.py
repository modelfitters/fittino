import os

path_to_slha1 = "SPheno.spc"
output_path = "SPheno.spc"
import models.check_param_card as converter
converter.convert_to_mg5card(path_to_slha1, output_path)
