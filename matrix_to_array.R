matrix_to_array = function(matrix,K, numSample, numTime){

  output = array(NA, dim = c(K, numTime, numSample))
  for (i in 1:numSample){
    output[,,i]= matrix[,(numTime*i-(numTime-1)):(numTime*i)]
  }
  return(output)
}
