from datetime import datetime
from typing import Optional, List
from pydantic import BaseModel


# BaseModel을 사용해야함...?
# 이게 engine 하고의 연결, 결국은 json 받어 받는데 어떻게 받음? BaseModel에 받는게 아니고
# 그냥 받어서 Model에 direct로 넣어 이게 제일 깔끔한데...? 중간에 그런 Adapter가 필요한가..?

class IRModel(BaseModel):
    name: str


class G2Model(IRModel):
    pass


class ModelCalibration(BaseModel):
    evaluationDate: datetime
    calendar: str
    vol_matrix: List[List]
    option_maturities: List[str]
    swap_maturity: List[str]
    model: IRModel
    # description: Optional[str] = None


model = ModelCalibration(
    evaluationDate=datetime(2021,4,30),
    calendar='kr'
    model=IRModel(name='testcalibration'))


print(model.dict())