# Estimation-of-RUL

This project was aimed at estimating the residual useful life of a Li-ion battery with the Constant Current Constant Voltage(CC-CV) charging mechanism. The datasets used in this project are inspired from the datasets available on Prognostic Centre of Excellence.

The project is phased out in two phases:

1. **First principle** approach where the working of battery is studied from the electrochemical prospective. In this approach, the battery is divided into finite elements and then electro chemical simulations are performed to estimate when the capacity of the battery reaches 80% its initial capacity

The results from the first principle approach validates the expected profile for the current, voltage and discharge capacity. 

<img width="415" alt="Screenshot 2023-05-21 215929" src="https://github.com/JavalVyas2000/Estimation-of-RUL/assets/73403218/ffea3c6c-4f15-4319-ad0a-d748dc9b1c25">


<img width="417" alt="image" src="https://github.com/JavalVyas2000/Estimation-of-RUL/assets/73403218/e8dfd5e6-3eeb-43e7-8020-902d286c32c4">


<img width="433" alt="image" src="https://github.com/JavalVyas2000/Estimation-of-RUL/assets/73403218/f50881e1-1bd9-45ea-a9c9-cdc29b8c2126">

2. **Data Driven** approach was used where the residual useful life of a battery was calculated based on the historical data recorded. This made the estimation to be fast in order of seconds and thus could be used for online purposes. 
This approach makes use of the following 3 methodologies:
a. **ARIMA model** 

<img width="431" alt="image" src="https://github.com/JavalVyas2000/Estimation-of-RUL/assets/73403218/2474c520-d727-4baa-aeed-02a756cef8cf">

b. **Fourier Series**

<img width="435" alt="image" src="https://github.com/JavalVyas2000/Estimation-of-RUL/assets/73403218/a643ea8f-6bc8-456b-95c6-4c98206b46d1">

c. **LSTM model**

<img width="434" alt="image" src="https://github.com/JavalVyas2000/Estimation-of-RUL/assets/73403218/33874d22-9670-4600-b37b-a50dc90ef5fc">


The results from the data driven methods are promising and specially from the ARIMA model. 
