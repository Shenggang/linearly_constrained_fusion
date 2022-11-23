consumption_raw = read.csv(file='consumption.csv')
survey_raw = read.csv(file='survey_useful_100.csv')

consumption_array = array(NaN, c(50, 535, 3))
survey_50 = array(1, c(50, 8))

count = 1
pointer = 1
while (count <= 50)
{
  if (survey_raw[pointer,2] > 0) # number of adults living at home
  {
    consumption_array[count,,] = array(data.matrix(consumption_raw[1:535+(pointer-1)*535, 4:6]), c(535,3))
    survey_50[count,] = array(data.matrix(survey_raw[pointer, 2:9]),c(1,8))
    count = count + 1
  }
  pointer = pointer+1
}

save(consumption_array, survey_50, file = "data_true.RData")