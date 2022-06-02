from datetime import datetime
from typing import Optional, List
import fastapi
from pydantic import BaseModel

# 기능?


# function 연결...? -> api를 두번부르는게 아닌 정보를 가지고 있다가 내부에서 
# 여러번 불러주는 형태...?
# 이거의 연결은 BaseModel 과 BaseModel을 연결함
# job의 경우는 기본적으로 BaseModel이 하나씩 associate되어 있어야하며
# 그 안에 class 는 model class 인 경우와 BaseModel을 조합해서 사용해야함


# instruments를 이용해 모델의 parameter를 구함
# 방법론
# 이놈은 조합되어서 들어옴. job 마다 BaseModel이 정의된거는 있어야겠네,
# 근데 그거의 item에 model class 를 넣을 수 있나...? 그냥 dict로 받아서 check 함.
# 있는 거는 
# job을 

# ql 결국은 여기로 박힘
# { tenors: ['3m'], indexes, calendar }

# output에 대한 정의?
class SwaptionSelect(BaseModel):
    vol_matrix: List[List]
    option_maturities: List[str]
    swap_maturity: List[str]

    # [{swap_maturity: '3m'}, {swap_maturity: '6m'}]
    def output(self):
        pass


 
# 여기는 ql이 필요한가...?
class ModelCalibration(BaseModel):
    evaluationDate: datetime
    calendar: str
    _swaptions: SwaptionSelect
    swaptions: List[str] # [ { swap_maturity: '3m'}, { swap_maturity: '6m'} ]
    model: IRModel



@fastapi
def model_calibration(ModelCalibration):

    pass


# market matrix를 받아서 picking함
# swaption instruments
@fastapi
def swaption_select():
    pass



# 조합...?
