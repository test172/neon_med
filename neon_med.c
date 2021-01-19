//*******************************************************************************************************
// Begin 添加的中值滤波的代码
//-------------------------------------------------------------------------------------------------------------------------
#define UINT16_MAX 65535


//-------------------------------------------------------------------------------------------------------------------------
/* Macros */

/* In-place swap top 64 bits of "a" with bottom 64 bits of "b" -- one operation (vswp)
*/
// 先把a和b的低四个数据和高四个数据分离开来（使用neon函数vget_low_u16和vget_high_u16）
// 然后把a的低四位数据和b的低四位数据结合起来（使用neon函数vcombine_u16）
// 把a的高四位数据和b的高四位数据结合起来
#define vtrn64q(a, b) \
  do { \
	uint16x4_t vtrn64_tmp_a0 = vget_low_u16(a), vtrn64_tmp_a1 = vget_high_u16(a); \
	uint16x4_t vtrn64_tmp_b0 = vget_low_u16(b), vtrn64_tmp_b1 = vget_high_u16(b); \
	{ \
	  uint16x4_t vtrn64_tmp = vtrn64_tmp_a1; \
	  vtrn64_tmp_a1 =  vtrn64_tmp_b0; \
	  vtrn64_tmp_b0 = vtrn64_tmp; \
	} \
	(a) = vcombine_u16(vtrn64_tmp_a0, vtrn64_tmp_a1); \
	(b) = vcombine_u16(vtrn64_tmp_b0, vtrn64_tmp_b1); \
     } while (0) 

/* In-place exchange odd 32-bit words of "a" with even 32-bit words of "b" -- one operation
*/
// 把a和b的16位的数据合成32位的数据（使用neon函数vreinterpretq_u32_u16）
// 使用vtrnq_u32得到组合好的a和b
// 重新把32为数据分解成16位的数据（使用neon函数vreinterpretq_u16_u32）
#define vtrn32q(a, b) \
  do  { \
	uint32x4x2_t vtrn32_tmp = vtrnq_u32(vreinterpretq_u32_u16(a), vreinterpretq_u32_u16(b)); \
	(a) = vreinterpretq_u16_u32(vtrn32_tmp.val[0]); \
	(b) = vreinterpretq_u16_u32(vtrn32_tmp.val[1]); \
      } while (0)

/* In-place exchange odd 16-bit words of "a" with even 16-bit words of "b" -- one operation 
*/
// 使用vtrnq_u16得到组合好的a和b
#define vtrn16q(a, b) \
  do { \
	uint16x8x2_t vtrn16_tmp = vtrnq_u16((a), (b)); \
	(a) = vtrn16_tmp.val[0]; \
	(b) = vtrn16_tmp.val[1]; \
     } while (0)

// 一开始： 排好序的数据是按照a的第一个数据<b的第一个数据<a的第二个数据<b的第二个数据<......
// 处理后： 使用vzipq_u16使得顺序为a的第一个数据<a的第二个数据<......<b的第一个数据<b的第二个数据<......
#define vzipq(a, b) \
  do { \
	uint16x8x2_t vzip_tmp = vzipq_u16((a), (b)); \
	(a) = vzip_tmp.val[0]; \
	(b) = vzip_tmp.val[1]; \
     } while (0)

// 两行数据(行a, 行b)进行比较，小的数据放在行a（使用neon函数vminq_u16）
// 两行数据(行a, 行b)进行比较，大的数据放在行b（使用neon函数vmaxq_u16）
#define vminmaxq(a, b) \
  do { \
	uint16x8_t minmax_tmp = (a); \
	(a) = vminq_u16((a), (b)); \
	(b) = vmaxq_u16(minmax_tmp, (b)); \
     } while (0)

// 使用vget_high_u16得到a的高的4个数据，然后使用vrev64_u16旋转180度
// 使用vget_low_u16得到a的低的4个数据，然后使用vrev64_u16旋转180度
// 把a的高四位数据和a的低四位数据结合起来（使用neon函数vcombine_u16）
#define vrev128q_u16(a) \
	vcombine_u16(vrev64_u16(vget_high_u16(a)), vrev64_u16(vget_low_u16(a)))

// 对一个7行8列的块，每一列排好了序
// spitch 是行的长度(640)
void loadblock(uint16_t dst[8][8], uint16_t const *src, int spitch)
{
// 5行排序
////////////////////////////////////////////////////////////////////////////////////////////////////////
#if 1
	// 5行数据（多取了3列）
	uint16x8_t q0, q1, q2, q3, q4, q5, q6, q7;

	// 从内存加载5行数据到neon寄存器
	//spitch >>=1;
	q0 = vld1q_u16(src); src += spitch;//指向下一行
	q1 = vld1q_u16(src); src += spitch;
	q2 = vld1q_u16(src); src += spitch;
	q3 = vld1q_u16(src); src += spitch;
	q4 = vld1q_u16(src); 

	/* vminmaxq() performs a pair of vminq and vmaxq operations to put the two arguments in order */
	// 对列进行排序，得到排好序的每一列数据

	// 相邻的两列排序，最后一行保持不变
	vminmaxq(q0, q1);
	vminmaxq(q2, q3);

	// 隔一列进行排序，第6行保持不变
	vminmaxq(q0, q2);
	vminmaxq(q1, q3);

	vminmaxq(q1, q2);
	
	vminmaxq(q3, q4);
	vminmaxq(q2, q3);
	vminmaxq(q1, q2);
	vminmaxq(q0, q1);



	// 生成全是最大值的3行(5行8列变成了8行8列)（多取了1列）
	q5 = vdupq_n_u16(UINT16_MAX);
	q6 = vdupq_n_u16(UINT16_MAX);

#endif
////////////////////////////////////////////////////////////////////////////////////////////////////////

	// 生成全是最大值的1行(5行8列变成了8行8列)（多取了1列）
	q7 = vdupq_n_u16(UINT16_MAX);
	
	
	// 
	// 处理后： 排的数据是按照a的第一个数据<b的第一个数据<a的第二个数据<b的第二个数据<......（使用定义的宏vzipq）
	vzipq(q0, q1);
	vzipq(q2, q3);
	vzipq(q4, q5);
	//vzipq(q6, q7);

	uint32x4x4_t tmp;
	// 把16位的数据合成32位的数据（使用neon函数vreinterpretq_u32_u16）
	tmp.val[0] = vreinterpretq_u32_u16(q0); tmp.val[1] = vreinterpretq_u32_u16(q2);
	tmp.val[2] = vreinterpretq_u32_u16(q4); tmp.val[3] = vreinterpretq_u32_u16(q6);

	// 把4个寄存器的数据放到内存（使用neon函数vst4q_u32）
	vst4q_u32((uint32_t *)&dst[0], tmp);

	// 把16位的数据合成32位的数据（使用neon函数vreinterpretq_u32_u16）
	tmp.val[0] = vreinterpretq_u32_u16(q1); tmp.val[1] = vreinterpretq_u32_u16(q3);
	tmp.val[2] = vreinterpretq_u32_u16(q5); tmp.val[3] = vreinterpretq_u32_u16(q7);

	// 把4个寄存器的数据放到内存（使用neon函数vst4q_u32）
	vst4q_u32((uint32_t *)&dst[4], tmp);	
}
//-------------------------------------------------------------------------------------------------------------------------


// 对两行数据进行排序
static inline uint16x8x2_t bitonic_resort_16(uint16x8_t a, uint16x8_t b)
{
	/* Start with two vectors:
	 * +---+---+---+---+---+---+---+---+
	 * | a | b | c | d | e | f | g | h |
	 * +---+---+---+---+---+---+---+---+
	 * +---+---+---+---+---+---+---+---+
	 * | i | j | k | l | m | m | o | p |
	 * +---+---+---+---+---+---+---+---+
	 * All the elements of the first are guaranteed to be less than or equal to
	 * all of the elements in the second, and both vectors are bitonic.
	 * We need to perform these operations to completely sort both lists:
	 * 行a的左边的四个和右边的四个进行比较,(结果)小的数据放在行a的左边四个，大的数据放在行b的左边四个
	 * 行b的左边的四个和右边的四个进行比较,(结果)小的数据放在行a的右边四个，大的数据放在行b的右边四个
	 *	vminmax([abcd],[efgh])		vminmax([ijkl],[mnop]) 
	 */

	/* We can align the necessary pairs for mixing with transpose operations, 
	 * like so:
	 */

	// 1.重新排列行a和行b
	vtrn64q(a,b);
	/* We now have:	
	* +---+---+---+---+---+---+---+---+
	* | a | b | c | d | i | j | k | l |
	* +---+---+---+---+---+---+---+---+
	* +---+---+---+---+---+---+---+---+
	* | e | f | g | h | m | m | o | p |
	* +---+---+---+---+---+---+---+---+
	*/
	// 两行数据(行a, 行b)进行比较，小的数据放在行a, 大的数据放在行b
	vminmaxq(a,b);
	// 2.重新排列行a和行b
	vtrn32q(a,b);
	/* We now have:	
	* +---+---+---+---+---+---+---+---+
	* | a | b | e | f | i | j | m | m |
	* +---+---+---+---+---+---+---+---+
	* +---+---+---+---+---+---+---+---+
	* | c | d | g | h | k | l | o | p |
	* +---+---+---+---+---+---+---+---+
	*/
	// 两行数据(行a, 行b)进行比较，小的数据放在行a, 大的数据放在行b
	vminmaxq(a,b);
	// 3.重新排列行a和行b
	vtrn16q(a,b);
	/* We now have:	
	* +---+---+---+---+---+---+---+---+
	* | a | c | e | g | i | k | m | o |
	* +---+---+---+---+---+---+---+---+
	* +---+---+---+---+---+---+---+---+
	* | b | d | f | h | j | l | n | p |
	* +---+---+---+---+---+---+---+---+
	*/
	// 两行数据(行a, 行b)进行比较，小的数据放在行a, 大的数据放在行b
	vminmaxq(a,b);
	/* Since we now have separate vectors of odd and even lanes, we can simply 
	* use vzip to stich them back together producing our original
	* positioning. Then we're done.
	*/
	// 把a和b按序连接起来
	return vzipq_u16(a,b);
}
//-------------------------------------------------------------------------------------------------------------------------




static inline uint16x8x2_t bitonic_merge_16(uint16x8_t a, uint16x8_t b)
{
	// 翻转行b(变成了从大到小排列)
	b = vrev128q_u16(b);
	// 两行数据(行a, 行b)进行比较，小的数据放在行a, 大的数据放在行b
	// 数据分成了两部分
	vminmaxq(a,b);
	// 排序
	return bitonic_resort_16(a,b);
}




//-------------------------------------------------------------------------------------------------------------------------
/* Macros */
// 两行数据(行a, 行b)进行比较，大的数据放在行a（使用neon函数vmaxq_u16）
// 两行数据(行a, 行b)进行比较，小的数据放在行b（使用neon函数vminq_u16）
#define vmaxminq(a, b) \
	do { \
	uint16x8_t maxmin_tmp = (a); \
	(a) = vmaxq_u16((a), (b)); \
	(b) = vminq_u16(maxmin_tmp, (b)); \
	} while (0)

// 对一行数据进行排序
#define vmaxmin_half(a) \
	do { \
	uint16x4_t minmax_tmp_lo = vget_low_u16(a), minmax_tmp_hi = vget_high_u16(a); \
	vmaxmin(minmax_tmp_lo, minmax_tmp_hi); \
	(a) = vcombine_u16(minmax_tmp_lo, minmax_tmp_hi); \
	} while (0)

// 大的四个数据放在前面四个，小的四个数据放在后面四个
#define vmaxmin(a, b) \
	do { \
	uint16x4_t maxmin_tmp = (a); \
	(a) = vmax_u16((a), (b)); \
	(b) = vmin_u16(maxmin_tmp, (b)); \
	} while (0)


///////////////////////////////////////////////////////////////////////////////////////
// 把a和b的16位的数据合成32位的数据（使用neon函数vreinterpretq_u32_u16）
// 使用vtrnq_u32得到组合好的a和b
// 重新把32为数据分解成16位的数据（使用neon函数vreinterpretq_u16_u32）
#define vtrn32(a, b) \
  do  { \
	uint32x2x2_t vtrn32_tmp = vtrn_u32(vreinterpret_u32_u16(a), vreinterpret_u32_u16(b)); \
	(a) = vreinterpret_u16_u32(vtrn32_tmp.val[0]); \
	(b) = vreinterpret_u16_u32(vtrn32_tmp.val[1]); \
      } while (0)

/* In-place exchange odd 16-bit words of "a" with even 16-bit words of "b" -- one operation 
*/
// 使用vtrn_u16得到组合好的a和b
// implement by myself
#define vtrn16(a, b) \
  do { \
	uint16x4x2_t vtrn16_tmp = vtrn_u16((a), (b)); \
	(a) = vtrn16_tmp.val[0]; \
	(b) = vtrn16_tmp.val[1]; \
     } while (0)

// 一开始： 排好序的数据是按照a的第一个数据<b的第一个数据<a的第二个数据<b的第二个数据<......
// 处理后： 使用vzipq_u16使得顺序为a的第一个数据<a的第二个数据<......<b的第一个数据<b的第二个数据<......
// implement by myself
#define vzip(a, b) \
  do { \
	uint16x4x2_t vzip_tmp = vzip_u16((a), (b)); \
	(a) = vzip_tmp.val[0]; \
	(b) = vzip_tmp.val[1]; \
     } while (0)
///////////////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////////////
// 对两行数据进行排序(根据bitonic_resort_16修改的)
static inline uint16x8x2_t bitonic_resort_16r(uint16x8_t a, uint16x8_t b)
{
	/* Start with two vectors:
	 * +---+---+---+---+---+---+---+---+
	 * | a | b | c | d | e | f | g | h |
	 * +---+---+---+---+---+---+---+---+
	 * +---+---+---+---+---+---+---+---+
	 * | i | j | k | l | m | m | o | p |
	 * +---+---+---+---+---+---+---+---+
	 * All the elements of the first are guaranteed to be less than or equal to
	 * all of the elements in the second, and both vectors are bitonic.
	 * We need to perform these operations to completely sort both lists:
	 * 行a的左边的四个和右边的四个进行比较,(结果)小的数据放在行a的左边四个，大的数据放在行b的左边四个
	 * 行b的左边的四个和右边的四个进行比较,(结果)小的数据放在行a的右边四个，大的数据放在行b的右边四个
	 *	vminmax([abcd],[efgh])		vminmax([ijkl],[mnop]) 
	 */

	/* We can align the necessary pairs for mixing with transpose operations, 
	 * like so:
	 */

	// 1.重新排列行a和行b
	vtrn64q(a,b);
	/* We now have:	
	* +---+---+---+---+---+---+---+---+
	* | a | b | c | d | i | j | k | l |
	* +---+---+---+---+---+---+---+---+
	* +---+---+---+---+---+---+---+---+
	* | e | f | g | h | m | m | o | p |
	* +---+---+---+---+---+---+---+---+
	*/
	// 两行数据(行a, 行b)进行比较，大的数据放在行a, 小的数据放在行b
	vmaxminq(a,b);
	// 2.重新排列行a和行b
	vtrn32q(a,b);
	/* We now have:	
	* +---+---+---+---+---+---+---+---+
	* | a | b | e | f | i | j | m | m |
	* +---+---+---+---+---+---+---+---+
	* +---+---+---+---+---+---+---+---+
	* | c | d | g | h | k | l | o | p |
	* +---+---+---+---+---+---+---+---+
	*/
	// 两行数据(行a, 行b)进行比较，大的数据放在行a, 小的数据放在行b
	vmaxminq(a,b);
	// 3.重新排列行a和行b
	vtrn16q(a,b);
	/* We now have:	
	* +---+---+---+---+---+---+---+---+
	* | a | c | e | g | i | k | m | o |
	* +---+---+---+---+---+---+---+---+
	* +---+---+---+---+---+---+---+---+
	* | b | d | f | h | j | l | n | p |
	* +---+---+---+---+---+---+---+---+
	*/
	// 两行数据(行a, 行b)进行比较，大的数据放在行a, 小的数据放在行b
	vmaxminq(a,b);
	/* Since we now have separate vectors of odd and even lanes, we can simply 
	* use vzip to stich them back together producing our original
	* positioning. Then we're done.
	*/
	// 把a和b按序连接起来
	return vzipq_u16(a,b);
}
///////////////////////////////////////////////////////////////////////////////////////

static inline uint16x8_t bitonic_resort_8r(uint16x8_t a)
{
	uint16x4_t a0 = vget_low_u16(a), a1 = vget_high_u16(a);
	
	vtrn32(a0, a1);
	vmaxmin(a0, a1);
	vtrn16(a0, a1);
	vmaxmin(a0, a1);
	vzip(a0, a1);
	return vcombine_u16(a0, a1);
}

// 排好序的三行数据，按从大到小的顺序排的(a是排好序的)
static inline uint16x8x3_t bitonic_merge_24r(uint16x8x2_t a, uint16x8_t b)
{
	uint16x8x3_t result;

	// 行b旋转180度（vrev128q_u16是定义的宏）,变成了从大到小排序
	b = vrev128q_u16(b);	
	
	//三行对应的最大值放在a.val[0], 最小值放到了b行中
	vmaxminq(a.val[0], b);
	vmaxminq(a.val[0], a.val[1]);

	// 得到排好序的两行?????????????????????????????????????? name may be wrong
	a = bitonic_resort_16r(a.val[0], a.val[1]);

	// 对行b进行排序
	vmaxmin_half(b);
	b = bitonic_resort_8r(b);

	result.val[0] = a.val[0];
	result.val[1] = a.val[1];
	result.val[2] = b;

	return result;
}
//如何merge到32-lane的向量？40-lane的数据如何取中值？
//例子显示是2个16-lane的数据merge到1个32-lane的数据。
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*?????????????????????????????????????????????????????????????????????????????????????????????????????????*/
// 总共是5行8列， 中值在第3行的第1个值(排好序的a,排好序的b0,排好序的b1)
// 这个实现有没有问题 （implement by myself）
// 中值是第2行的第5个
static uint16x8x2_t bitonic_median_40(uint16_t *dst, uint16x8x2_t a, uint16x8x2_t b0, uint16x8_t b1)
{
	// merge 3组vector(完全排好序的三组vector)按从大到小排列
	// 第一列全是65535
	uint16x8x3_t b = bitonic_merge_24r(b0, b1); /* MUST inline */
	

	b.val[1] = vminq_u16(a.val[0], b.val[1]);
	b.val[2] = vmaxq_u16(b.val[1], b.val[2]);
        //TO DO
        //还要和a的最后一行比较下得到最小值
	{
	  uint16x8_t tmp;
	  vmaxmin_half(b.val[2]);
	  tmp = bitonic_resort_8r(b.val[2]);
	  // dst是中值
	  vst1q_lane_u16(dst, tmp, 3);
	}
	return a;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////




//-------------------------------------------------------------------------------------------------------------------------
	

// count是图像的总像素（640*480）
// spitch 是行的长度(640)
void filter_row_bs(uint16_t *dst, uint16_t const *src, int spitch, int count)
{
	dst += 640*2+2;
	uint16_t scratch[8][8];
	struct
	{
	  uint16x8x2_t ef;
	  uint16x8_t b, d, h;
	} state;
	// scratch 是每一行都排好序的向量
	// 8行数据(多排了3行) src是8列8列的取，从左往右，每次都多取了3列
	loadblock(scratch, src, spitch);
	
	
	uint16x8_t c, d;
	//把R通道的值提取出来
	c = vld1q_u16(scratch[1]);
	// 第三行
	d = vld1q_u16(scratch[2]);

	// 第二，三行合在一起
	state.ef = bitonic_merge_16(c, d);
	  
	// 第一行
	state.b = vld1q_u16(scratch[0]);
	// 倒数第二行
	state.h = vld1q_u16(scratch[3]);
	// 第三行
	state.d = d;
	// 指向图像第5列
	src += 4;

	/* We maintain three sets (components) of partially combined lists and we recombine
	*  these for a full bitonic sort in various permutations in an eight-phase rota.
	*
	* Register allocation collapses in a heap, here, and we spend a lot of our
	* time marshalling things through the stack. We could probably do much
	* better with hand-written assembly.
	*/

	while (count > 0)
	{
	  //int i;
	  //uint16_t const *pft, *pft2;

	  // scratch 是每一行都排好序的向量, 总共有24行
	  // 再次加载src到scratch
	  // 这里是RGB三个通道的值
	  // 对于R通道就是8行值
	  loadblock(scratch, src, spitch);
	  src += 8;
	  // 指向图像第5列+8列
	  //pft = src;
	  // 最后列数少于8列
	  //if (count < 8) pft -= 8;
	  //pft2 = pft+4*spitch;//指向图像第5行

	  
	  //int cnt = count;
	  uint16x8x2_t ef;
	  uint16x8_t h;

	  // x y轮流得到最后一行
	  uint16x8_t x, y;
	  uint16x8x2_t xx;	

	  // 第一个数据块的ef
	  ef = state.ef;
		
		
	  /* stage 0 */
	  x = vld1q_u16(scratch[0]);
	  //__builtin_prefetch(pft+0);
	  // 最后两行，会重新使用三次
	  
	  xx = bitonic_merge_16(state.h, x);

	  //uint16_t ttt[8];
	  
	  // 作为下一次重新使用（行5）
	  h = x;
	  // 第一个中值，放在dst里
	  // ef还要重新使用一次，xx还要重新使用三次,state.b使用完了
	  xx = bitonic_median_40(dst, xx, ef, state.b);
	  //printf("%d\n", *dst);
	  //if (--cnt <= 0)
		//continue;
	  #if 1
	  /* stage 1 */
	  y = vld1q_u16(scratch[1]);
	  //__builtin_prefetch(pft+16);
	  // ef重新使用完了，xx还要重新使用两次
	  xx = bitonic_median_40(dst+1, xx, ef, y);
	  //if (--cnt <= 0)
		//continue;

	  /* stage 2 */
	  x = vld1q_u16(scratch[2]);
	  //__builtin_prefetch(pft+23); pft += spitch;
	  // 最后两行，会重新使用三次
	  
	  ef = bitonic_merge_16(y, x);
	  // 作为下一次重新使用（行7）
	  state.b = x;
	  //ef = xx;
	  // ef重新使用三次，xx还要重新使用一次,state.d使用完了
	  ef = bitonic_median_40(dst+2, ef, xx, state.d);
	  //if (--cnt <= 0)
		//continue;

	  /* stage 3 */
	  y = vld1q_u16(scratch[3]);
	  //__builtin_prefetch(pft+0);
	  // xx重新使用完了，ef还要重新使用两次
	  
	  ef = bitonic_median_40(dst+3, ef, xx, y);
	  //if (--cnt <= 0)
		//continue;

	  /* stage 4 */
	  x = vld1q_u16(scratch[4]);
	  
	  //__builtin_prefetch(pft+16);
	  // 最后两行，会重新使用三次
	  xx = bitonic_merge_16(y, x);	  
	  // 作为下一次重新使用(行9)
	  state.d = x;
	  // ef还要重新使用一次，xx还要重新使用三次, h使用完了
	  xx = bitonic_median_40(dst+4, xx, ef, h);
	  //if (--cnt <= 0)
		//continue;

	  /* stage 5 */
	  y = vld1q_u16(scratch[5]);
	  //__builtin_prefetch(pft+23); pft += spitch;
	  // ef重新使用完了，xx还要重新使用两次
	  xx = bitonic_median_40(dst+5, xx, ef, y);
	  //if (--cnt <= 0)
	  	//continue;

	  /* stage 6 */
	  x = vld1q_u16(scratch[6]);
	  //__builtin_prefetch(pft2); pft2 += 8;
	  // ef重新使用三次，xx还要重新使用一次
	  ef = bitonic_merge_16(y, x);
	  // 保存到第一个数据块，留待下一次循环使用------------------------------------
	  state.ef = ef;
		
	  // ef重新使用三次，xx还要重新使用一次, state.b使用完了
	  ef = bitonic_median_40(dst+6, ef, xx, state.b);

	  state.b = state.d;
	  // 作为下一次重新使用(行11)
	  state.d = x;
	  //if (--cnt <= 0)
		//continue;
	
	  /* stage 7 */
	  y = vld1q_u16(scratch[7]);
	  state.h = y;
	  // ef重新使用两次，xx使用完了
	  (void)bitonic_median_40(dst+7, ef, xx, y);
	

	  // 中值
	  dst += 8;

	  // 像素总个数
	  count -= 8;
	}
	#endif
}
//*******************************************************************************************************
// End 添加的中值滤波的代码
