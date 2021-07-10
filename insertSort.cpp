void insertionSort(vector<int>& a){
    for(int i = 1; i < a.size(); i++){
        int v = a[i];
        int j = i-1;
        while(j >= 0 && a[j] > v){
            a[j+1] = a[j];
            j--;
        }
        a[j+1] = v;
    }
}
