template<typename T>
void selectionSort(vector<T>& a){
    int minj;
    for(int i = 0; i < a.size() - 1; i++){
        minj = i;
        for(int j = i; j < a.size(); j++){
            if(a[j] < a[minj]) minj = j;
        }
        swap(a[i], a[minj]);
    }
}
