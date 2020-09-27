public class GeoSum {

    public static void main(String[] args) {

    	
    	
        // 5 to power of 5 is 3125
        System.out.println(GeoSum(5, 5));

       
    }

    public static int GeoSum(int a, int N) {

        int answer = 1;

        for (int i = 0; i < N; i++) {


            answer = answer * N;
            
            

        }

        return answer;

    }

}