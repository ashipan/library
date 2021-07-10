template<typename T>
void bubbleSort(vector<T>& a){
    bool flag = 1;
    for(int i = 0; flag; i++){
        flag = 0;
        for(int j = a.size()-1; j >= i+1; j--){
            if(a[j] < a[j-1]){
                swap(a[j], a[j-1]);
                flag = 1;
            }
        }
    }
}
