package MyJama.util;

public class Maths {

   /** sqrt(a^2 + b^2) without under/overflow. **/

   public static float hypot(float a, float b) {
      float r;
      if (Math.abs(a) > Math.abs(b)) {
         r = b/a;
         r = Math.abs(a)*(float)Math.sqrt(1+r*r);
      } else if (b != 0) {
         r = a/b;
         r = Math.abs(b)*(float)Math.sqrt(1+r*r);
      } else {
         r = 0.0f;
      }
      return r;
   }
}
