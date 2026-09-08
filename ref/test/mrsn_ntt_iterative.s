	.file	"mrsn_ntt_iterative.c"
	.text
	.p2align 4
	.globl	mrsn_ntt_256
	.type	mrsn_ntt_256, @function
mrsn_ntt_256:
.LFB49:
	.cfi_startproc
	endbr64
	pushq	%r15
	.cfi_def_cfa_offset 16
	.cfi_offset 15, -16
	movq	%rdi, %r15
	pushq	%r14
	.cfi_def_cfa_offset 24
	.cfi_offset 14, -24
	pushq	%r13
	.cfi_def_cfa_offset 32
	.cfi_offset 13, -32
	pushq	%r12
	.cfi_def_cfa_offset 40
	.cfi_offset 12, -40
	pushq	%rbp
	.cfi_def_cfa_offset 48
	.cfi_offset 6, -48
	pushq	%rbx
	.cfi_def_cfa_offset 56
	.cfi_offset 3, -56
	subq	$48, %rsp
	.cfi_def_cfa_offset 104
	movl	256(%r15), %ecx
	movl	(%rdi), %edi
	movl	768(%r15), %edx
	movl	512(%r15), %esi
	movl	%ecx, %eax
	shrl	$16, %ecx
	sall	$15, %eax
	andl	$2147450880, %eax
	orl	%ecx, %eax
	movl	%edx, %ecx
	shrl	$16, %edx
	sall	$15, %ecx
	andl	$2147450880, %ecx
	orl	%edx, %ecx
	movl	%ecx, %r8d
	addl	%eax, %ecx
	xorl	$2147483647, %r8d
	addl	%eax, %r8d
	movl	%ecx, %eax
	shrl	$31, %ecx
	movl	%r8d, %edx
	shrl	$31, %r8d
	andl	$2147483647, %eax
	andl	$2147483647, %edx
	addl	%ecx, %eax
	addl	%r8d, %edx
	leal	(%rdi,%rdx), %ecx
	notl	%edx
	movl	%ecx, %r8d
	andl	$2147483647, %ecx
	andl	$2147483647, %edx
	shrl	$31, %r8d
	addl	%edi, %edx
	addl	%r8d, %ecx
	movl	%ecx, (%r15)
	leal	(%rsi,%rax), %ecx
	notl	%eax
	movl	%ecx, %r8d
	shrl	$31, %ecx
	andl	$2147483647, %eax
	andl	$2147483647, %r8d
	addl	%esi, %eax
	addl	%r8d, %ecx
	leaq	omegas512(%rip), %r8
	movl	%ecx, 256(%r15)
	movl	%edx, %ecx
	andl	$2147483647, %edx
	movq	%r8, %rdi
	shrl	$31, %ecx
	leaq	496(%r8), %rsi
	addl	%ecx, %edx
	leaq	4(%r15), %rcx
	movl	%edx, 512(%r15)
	movl	%eax, %edx
	andl	$2147483647, %edx
	shrl	$31, %eax
	addl	%edx, %eax
	movl	%eax, 768(%r15)
	.p2align 4,,10
	.p2align 3
.L2:
	movl	(%rcx), %ebx
	movl	(%rdi), %r11d
	movl	4(%rdi), %r9d
	movl	256(%rcx), %r13d
	movq	%rbx, %r10
	movl	4(%rsi), %r12d
	imulq	%r11, %r10
	movq	%r10, %rax
	andl	$2147483647, %r10d
	shrq	$31, %rax
	addl	%eax, %r10d
	movl	512(%rcx), %eax
	movq	%rax, %rdx
	imulq	%r11, %rax
	imulq	%r9, %rdx
	imulq	%rbx, %r9
	movq	%rax, %r11
	andl	$2147483647, %eax
	movq	%rdx, %rbp
	andl	$2147483647, %edx
	shrq	$31, %r11
	shrq	$31, %rbp
	addl	%r11d, %eax
	movq	%r9, %rbx
	andl	$2147483647, %r9d
	addl	%ebp, %edx
	shrq	$31, %rbx
	movl	%edx, %r11d
	andl	$2147483647, %edx
	addl	%ebx, %r9d
	shrl	$31, %r11d
	addl	%r11d, %edx
	movl	%r10d, %r11d
	andl	$2147483647, %r10d
	shrl	$31, %r11d
	notl	%edx
	addl	%r11d, %r10d
	andl	$2147483647, %edx
	movl	(%rsi), %r11d
	addl	%r10d, %edx
	movl	%r9d, %r10d
	shrl	$31, %r9d
	movl	%edx, %ebx
	andl	$2147483647, %edx
	andl	$2147483647, %r10d
	shrl	$31, %ebx
	addl	%edx, %ebx
	movl	%eax, %edx
	andl	$2147483647, %edx
	addl	%r10d, %edx
	addl	%r9d, %edx
	shrl	$31, %eax
	movq	%rsi, %r9
	addl	%edx, %eax
	movq	%r13, %rdx
	imulq	%r12, %rdx
	movl	%eax, %r10d
	andl	$2147483647, %eax
	shrl	$31, %r10d
	addl	%eax, %r10d
	movq	%rdx, %rax
	andl	$2147483647, %edx
	shrq	$31, %rax
	addl	%eax, %edx
	movl	768(%rcx), %eax
	movq	%rax, %rbp
	imulq	%r12, %rax
	imulq	%r11, %rbp
	imulq	%r13, %r11
	movq	%rax, %r12
	andl	$2147483647, %eax
	movq	%rbp, %r14
	andl	$2147483647, %ebp
	shrq	$31, %r12
	shrq	$31, %r14
	addl	%r12d, %eax
	movq	%r11, %r13
	andl	$2147483647, %r11d
	addl	%r14d, %ebp
	shrq	$31, %r13
	movl	%ebp, %r12d
	andl	$2147483647, %ebp
	addl	%r13d, %r11d
	shrl	$31, %r12d
	addl	%r12d, %ebp
	movl	%edx, %r12d
	andl	$2147483647, %edx
	shrl	$31, %r12d
	notl	%ebp
	addl	%r12d, %edx
	andl	$2147483647, %ebp
	movl	%r11d, %r12d
	addl	%edx, %ebp
	movl	%ebp, %edx
	andl	$2147483647, %ebp
	shrl	$31, %edx
	addl	%ebp, %edx
	movl	%eax, %ebp
	andl	$2147483647, %ebp
	andl	$2147483647, %r12d
	shrl	$31, %r11d
	addq	$4, %rcx
	addl	%r12d, %ebp
	shrl	$31, %eax
	addq	$8, %rdi
	subq	$8, %rsi
	addl	%ebp, %r11d
	addl	%eax, %r11d
	movl	%r11d, %eax
	shrl	$31, %r11d
	andl	$2147483647, %eax
	addl	%r11d, %eax
	leal	(%rdx,%rbx), %r11d
	notl	%edx
	movl	%r11d, %ebp
	shrl	$31, %r11d
	andl	$2147483647, %edx
	andl	$2147483647, %ebp
	addl	%ebx, %edx
	addl	%ebp, %r11d
	movl	%r11d, -4(%rcx)
	leal	(%rax,%r10), %r11d
	notl	%eax
	movl	%r11d, %ebp
	shrl	$31, %r11d
	andl	$2147483647, %eax
	andl	$2147483647, %ebp
	addl	%r10d, %eax
	addl	%ebp, %r11d
	movl	%r11d, 252(%rcx)
	movl	%edx, %r11d
	shrl	$31, %edx
	andl	$2147483647, %r11d
	addl	%r11d, %edx
	movl	%edx, 508(%rcx)
	movl	%eax, %edx
	andl	$2147483647, %eax
	shrl	$31, %edx
	addl	%edx, %eax
	movl	%eax, 764(%rcx)
	cmpq	%r9, %r8
	jne	.L2
	movdqu	(%r15), %xmm5
	leaq	16(%r15), %rax
	movq	%r15, %r12
	movdqu	128(%r15), %xmm1
	movdqa	.LC0(%rip), %xmm0
	movq	%rax, 24(%rsp)
	movdqu	256(%r15), %xmm4
	movdqa	%xmm5, %xmm2
	paddd	%xmm1, %xmm2
	pandn	%xmm0, %xmm1
	movdqa	%xmm2, %xmm3
	pand	%xmm0, %xmm2
	paddd	%xmm5, %xmm1
	movdqu	16(%r15), %xmm5
	psrld	$31, %xmm3
	paddd	%xmm3, %xmm2
	movups	%xmm2, (%r15)
	movdqu	384(%r15), %xmm2
	paddd	%xmm4, %xmm2
	movdqa	%xmm2, %xmm3
	pand	%xmm0, %xmm2
	psrld	$31, %xmm3
	paddd	%xmm3, %xmm2
	movups	%xmm2, 128(%r15)
	movdqa	%xmm1, %xmm2
	pand	%xmm0, %xmm1
	psrld	$31, %xmm2
	paddd	%xmm2, %xmm1
	movups	%xmm1, 256(%r15)
	movdqu	384(%r15), %xmm1
	pandn	%xmm0, %xmm1
	paddd	%xmm4, %xmm1
	movdqu	272(%r15), %xmm4
	movdqa	%xmm1, %xmm2
	psrld	$31, %xmm1
	pand	%xmm0, %xmm2
	paddd	%xmm2, %xmm1
	movups	%xmm1, 384(%r15)
	movdqu	144(%r15), %xmm1
	movdqa	%xmm1, %xmm2
	pandn	%xmm0, %xmm1
	paddd	%xmm5, %xmm2
	paddd	%xmm5, %xmm1
	movdqu	32(%r15), %xmm5
	movdqa	%xmm2, %xmm3
	pand	%xmm0, %xmm2
	psrld	$31, %xmm3
	paddd	%xmm3, %xmm2
	movups	%xmm2, 16(%r15)
	movdqu	400(%r15), %xmm2
	paddd	%xmm4, %xmm2
	movdqa	%xmm2, %xmm3
	pand	%xmm0, %xmm2
	psrld	$31, %xmm3
	paddd	%xmm3, %xmm2
	movups	%xmm2, 144(%r15)
	movdqa	%xmm1, %xmm2
	psrld	$31, %xmm1
	pand	%xmm0, %xmm2
	paddd	%xmm2, %xmm1
	movups	%xmm1, 272(%r15)
	movdqu	400(%r15), %xmm1
	pandn	%xmm0, %xmm1
	paddd	%xmm4, %xmm1
	movdqu	288(%r15), %xmm4
	movdqa	%xmm1, %xmm2
	psrld	$31, %xmm1
	pand	%xmm0, %xmm2
	paddd	%xmm2, %xmm1
	movdqa	%xmm5, %xmm2
	movups	%xmm1, 400(%r15)
	movdqu	160(%r15), %xmm1
	paddd	%xmm1, %xmm2
	pandn	%xmm0, %xmm1
	movdqa	%xmm2, %xmm3
	pand	%xmm0, %xmm2
	paddd	%xmm5, %xmm1
	movdqu	48(%r15), %xmm5
	psrld	$31, %xmm3
	paddd	%xmm3, %xmm2
	movups	%xmm2, 32(%r15)
	movdqu	416(%r15), %xmm2
	paddd	%xmm4, %xmm2
	movdqa	%xmm2, %xmm3
	pand	%xmm0, %xmm2
	psrld	$31, %xmm3
	paddd	%xmm3, %xmm2
	movups	%xmm2, 160(%r15)
	movdqa	%xmm1, %xmm2
	pand	%xmm0, %xmm1
	psrld	$31, %xmm2
	paddd	%xmm2, %xmm1
	movups	%xmm1, 288(%r15)
	movdqu	416(%r15), %xmm1
	pandn	%xmm0, %xmm1
	paddd	%xmm4, %xmm1
	movdqu	304(%r15), %xmm4
	movdqa	%xmm1, %xmm2
	pand	%xmm0, %xmm1
	psrld	$31, %xmm2
	paddd	%xmm2, %xmm1
	movdqa	%xmm5, %xmm2
	movups	%xmm1, 416(%r15)
	movdqu	176(%r15), %xmm1
	paddd	%xmm1, %xmm2
	pandn	%xmm0, %xmm1
	movdqa	%xmm2, %xmm3
	psrld	$31, %xmm2
	paddd	%xmm5, %xmm1
	pand	%xmm0, %xmm3
	paddd	%xmm3, %xmm2
	movups	%xmm2, 48(%r15)
	movdqu	64(%r15), %xmm5
	movdqu	432(%r15), %xmm2
	paddd	%xmm4, %xmm2
	movdqa	%xmm2, %xmm3
	pand	%xmm0, %xmm2
	psrld	$31, %xmm3
	paddd	%xmm3, %xmm2
	movups	%xmm2, 176(%r15)
	movdqa	%xmm1, %xmm2
	pand	%xmm0, %xmm1
	psrld	$31, %xmm2
	paddd	%xmm2, %xmm1
	movups	%xmm1, 304(%r15)
	movdqu	432(%r15), %xmm1
	pandn	%xmm0, %xmm1
	paddd	%xmm4, %xmm1
	movdqu	320(%r15), %xmm4
	movdqa	%xmm1, %xmm2
	pand	%xmm0, %xmm1
	psrld	$31, %xmm2
	paddd	%xmm2, %xmm1
	movdqa	%xmm5, %xmm2
	movups	%xmm1, 432(%r15)
	movdqu	192(%r15), %xmm1
	paddd	%xmm1, %xmm2
	pandn	%xmm0, %xmm1
	movdqa	%xmm2, %xmm3
	pand	%xmm0, %xmm2
	paddd	%xmm5, %xmm1
	movdqu	80(%r15), %xmm5
	psrld	$31, %xmm3
	paddd	%xmm3, %xmm2
	movups	%xmm2, 64(%r15)
	movdqu	448(%r15), %xmm2
	paddd	%xmm4, %xmm2
	movdqa	%xmm2, %xmm3
	pand	%xmm0, %xmm2
	psrld	$31, %xmm3
	paddd	%xmm3, %xmm2
	movups	%xmm2, 192(%r15)
	movdqa	%xmm1, %xmm2
	pand	%xmm0, %xmm1
	psrld	$31, %xmm2
	paddd	%xmm2, %xmm1
	movups	%xmm1, 320(%r15)
	movdqu	448(%r15), %xmm1
	pandn	%xmm0, %xmm1
	paddd	%xmm4, %xmm1
	movdqu	336(%r15), %xmm4
	movdqa	%xmm1, %xmm2
	pand	%xmm0, %xmm1
	psrld	$31, %xmm2
	paddd	%xmm2, %xmm1
	movups	%xmm1, 448(%r15)
	movdqu	208(%r15), %xmm1
	movdqa	%xmm1, %xmm2
	pandn	%xmm0, %xmm1
	paddd	%xmm5, %xmm2
	paddd	%xmm5, %xmm1
	movdqu	96(%r15), %xmm5
	movdqa	%xmm2, %xmm3
	pand	%xmm0, %xmm2
	psrld	$31, %xmm3
	paddd	%xmm3, %xmm2
	movups	%xmm2, 80(%r15)
	movdqu	464(%r15), %xmm2
	paddd	%xmm4, %xmm2
	movdqa	%xmm2, %xmm3
	psrld	$31, %xmm2
	pand	%xmm0, %xmm3
	paddd	%xmm3, %xmm2
	movups	%xmm2, 208(%r15)
	movdqa	%xmm1, %xmm2
	psrld	$31, %xmm1
	pand	%xmm0, %xmm2
	paddd	%xmm2, %xmm1
	movups	%xmm1, 336(%r15)
	movdqu	464(%r15), %xmm1
	pandn	%xmm0, %xmm1
	paddd	%xmm4, %xmm1
	movdqu	352(%r15), %xmm4
	movdqa	%xmm1, %xmm2
	pand	%xmm0, %xmm1
	psrld	$31, %xmm2
	paddd	%xmm2, %xmm1
	movdqa	%xmm5, %xmm2
	movups	%xmm1, 464(%r15)
	movdqu	224(%r15), %xmm1
	paddd	%xmm1, %xmm2
	pandn	%xmm0, %xmm1
	movdqa	%xmm2, %xmm3
	psrld	$31, %xmm2
	paddd	%xmm5, %xmm1
	pand	%xmm0, %xmm3
	paddd	%xmm3, %xmm2
	movups	%xmm2, 96(%r15)
	movdqu	480(%r15), %xmm2
	paddd	%xmm4, %xmm2
	movdqa	%xmm2, %xmm3
	psrld	$31, %xmm2
	pand	%xmm0, %xmm3
	paddd	%xmm3, %xmm2
	movups	%xmm2, 224(%r15)
	movdqa	%xmm1, %xmm2
	psrld	$31, %xmm1
	pand	%xmm0, %xmm2
	paddd	%xmm2, %xmm1
	movups	%xmm1, 352(%r15)
	movdqu	480(%r15), %xmm1
	pandn	%xmm0, %xmm1
	paddd	%xmm4, %xmm1
	movdqa	%xmm1, %xmm2
	psrld	$31, %xmm1
	pand	%xmm0, %xmm2
	paddd	%xmm2, %xmm1
	movups	%xmm1, 480(%r15)
	movdqu	112(%r15), %xmm5
	movdqu	240(%r15), %xmm1
	movdqu	368(%r15), %xmm4
	movdqu	896(%r15), %xmm7
	movdqa	%xmm5, %xmm2
	paddd	%xmm1, %xmm2
	pandn	%xmm0, %xmm1
	movdqa	%xmm2, %xmm3
	pand	%xmm0, %xmm2
	paddd	%xmm5, %xmm1
	movdqu	768(%r15), %xmm5
	psrld	$31, %xmm3
	paddd	%xmm3, %xmm2
	movups	%xmm2, 112(%r15)
	movdqu	496(%r15), %xmm2
	paddd	%xmm4, %xmm2
	movdqa	%xmm2, %xmm3
	pand	%xmm0, %xmm2
	psrld	$31, %xmm3
	paddd	%xmm3, %xmm2
	movdqu	896(%r15), %xmm3
	movups	%xmm2, 240(%r15)
	movdqa	%xmm1, %xmm2
	pand	%xmm0, %xmm1
	psrld	$31, %xmm2
	pandn	%xmm0, %xmm3
	paddd	%xmm2, %xmm1
	movups	%xmm1, 368(%r15)
	movdqu	496(%r15), %xmm1
	pandn	%xmm0, %xmm1
	paddd	%xmm4, %xmm1
	movdqa	%xmm1, %xmm2
	pand	%xmm0, %xmm1
	psrld	$31, %xmm2
	paddd	%xmm2, %xmm1
	movdqu	512(%r15), %xmm2
	movups	%xmm1, 496(%r15)
	movdqu	640(%r15), %xmm1
	paddd	%xmm2, %xmm3
	paddd	%xmm7, %xmm2
	movdqa	%xmm3, %xmm4
	pand	%xmm0, %xmm3
	psrld	$31, %xmm4
	paddd	%xmm4, %xmm3
	movups	%xmm3, 512(%r15)
	movdqa	%xmm1, %xmm3
	pandn	%xmm0, %xmm1
	paddd	%xmm5, %xmm3
	paddd	%xmm5, %xmm1
	movdqu	528(%r15), %xmm5
	movdqa	%xmm3, %xmm4
	pand	%xmm0, %xmm3
	psrld	$31, %xmm4
	paddd	%xmm4, %xmm3
	movdqu	784(%r15), %xmm4
	movups	%xmm3, 640(%r15)
	movdqa	%xmm2, %xmm3
	pand	%xmm0, %xmm2
	psrld	$31, %xmm3
	paddd	%xmm3, %xmm2
	movups	%xmm2, 768(%r15)
	movdqa	%xmm1, %xmm2
	pand	%xmm0, %xmm1
	psrld	$31, %xmm2
	paddd	%xmm2, %xmm1
	movdqu	912(%r15), %xmm2
	movups	%xmm1, 896(%r15)
	movdqu	656(%r15), %xmm1
	pandn	%xmm0, %xmm2
	paddd	%xmm5, %xmm2
	movdqa	%xmm2, %xmm3
	pand	%xmm0, %xmm2
	psrld	$31, %xmm3
	paddd	%xmm3, %xmm2
	movups	%xmm2, 528(%r15)
	movdqa	%xmm1, %xmm2
	pandn	%xmm0, %xmm1
	paddd	%xmm4, %xmm2
	paddd	%xmm4, %xmm1
	movdqa	%xmm2, %xmm3
	pand	%xmm0, %xmm2
	psrld	$31, %xmm3
	paddd	%xmm3, %xmm2
	movups	%xmm2, 656(%r15)
	movdqu	912(%r15), %xmm2
	paddd	%xmm5, %xmm2
	movdqu	800(%r15), %xmm5
	movdqa	%xmm2, %xmm3
	pand	%xmm0, %xmm2
	psrld	$31, %xmm3
	paddd	%xmm3, %xmm2
	movdqu	928(%r15), %xmm3
	movups	%xmm2, 784(%r15)
	movdqa	%xmm1, %xmm2
	pand	%xmm0, %xmm1
	psrld	$31, %xmm2
	pandn	%xmm0, %xmm3
	paddd	%xmm2, %xmm1
	movdqu	544(%r15), %xmm2
	movups	%xmm1, 912(%r15)
	movdqu	672(%r15), %xmm1
	paddd	%xmm2, %xmm3
	movdqa	%xmm3, %xmm4
	pand	%xmm0, %xmm3
	psrld	$31, %xmm4
	paddd	%xmm4, %xmm3
	movups	%xmm3, 544(%r15)
	movdqa	%xmm1, %xmm3
	pandn	%xmm0, %xmm1
	paddd	%xmm5, %xmm3
	paddd	%xmm5, %xmm1
	movdqa	%xmm3, %xmm4
	pand	%xmm0, %xmm3
	psrld	$31, %xmm4
	paddd	%xmm4, %xmm3
	movups	%xmm3, 672(%r15)
	movdqu	928(%r15), %xmm7
	movdqu	816(%r15), %xmm5
	paddd	%xmm7, %xmm2
	movdqu	944(%r15), %xmm7
	movdqa	%xmm2, %xmm3
	pand	%xmm0, %xmm2
	psrld	$31, %xmm3
	paddd	%xmm3, %xmm2
	movdqu	944(%r15), %xmm3
	movups	%xmm2, 800(%r15)
	movdqa	%xmm1, %xmm2
	pand	%xmm0, %xmm1
	psrld	$31, %xmm2
	pandn	%xmm0, %xmm3
	paddd	%xmm2, %xmm1
	movdqu	560(%r15), %xmm2
	movups	%xmm1, 928(%r15)
	movdqu	688(%r15), %xmm1
	paddd	%xmm2, %xmm3
	paddd	%xmm7, %xmm2
	movdqu	960(%r15), %xmm7
	movdqa	%xmm3, %xmm4
	pand	%xmm0, %xmm3
	psrld	$31, %xmm4
	paddd	%xmm4, %xmm3
	movups	%xmm3, 560(%r15)
	movdqa	%xmm1, %xmm3
	pandn	%xmm0, %xmm1
	paddd	%xmm5, %xmm3
	paddd	%xmm5, %xmm1
	movdqu	832(%r15), %xmm5
	movdqa	%xmm3, %xmm4
	pand	%xmm0, %xmm3
	psrld	$31, %xmm4
	paddd	%xmm4, %xmm3
	movups	%xmm3, 688(%r15)
	movdqa	%xmm2, %xmm3
	pand	%xmm0, %xmm2
	psrld	$31, %xmm3
	paddd	%xmm3, %xmm2
	movdqu	960(%r15), %xmm3
	movups	%xmm2, 816(%r15)
	movdqa	%xmm1, %xmm2
	psrld	$31, %xmm1
	pand	%xmm0, %xmm2
	pandn	%xmm0, %xmm3
	paddd	%xmm2, %xmm1
	movdqu	576(%r15), %xmm2
	movups	%xmm1, 944(%r15)
	movdqu	704(%r15), %xmm1
	paddd	%xmm2, %xmm3
	paddd	%xmm7, %xmm2
	movdqa	%xmm3, %xmm4
	pand	%xmm0, %xmm3
	psrld	$31, %xmm4
	paddd	%xmm4, %xmm3
	movups	%xmm3, 576(%r15)
	movdqa	%xmm1, %xmm3
	pandn	%xmm0, %xmm1
	paddd	%xmm5, %xmm3
	paddd	%xmm5, %xmm1
	movdqu	592(%r15), %xmm5
	movdqa	%xmm3, %xmm4
	pand	%xmm0, %xmm3
	psrld	$31, %xmm4
	paddd	%xmm4, %xmm3
	movdqu	848(%r15), %xmm4
	movups	%xmm3, 704(%r15)
	movdqa	%xmm2, %xmm3
	pand	%xmm0, %xmm2
	psrld	$31, %xmm3
	paddd	%xmm3, %xmm2
	movups	%xmm2, 832(%r15)
	movdqa	%xmm1, %xmm2
	psrld	$31, %xmm1
	pand	%xmm0, %xmm2
	paddd	%xmm2, %xmm1
	movdqu	976(%r15), %xmm2
	movups	%xmm1, 960(%r15)
	movdqu	720(%r15), %xmm1
	pandn	%xmm0, %xmm2
	paddd	%xmm5, %xmm2
	movdqa	%xmm2, %xmm3
	psrld	$31, %xmm2
	pand	%xmm0, %xmm3
	paddd	%xmm3, %xmm2
	movups	%xmm2, 592(%r15)
	movdqa	%xmm4, %xmm2
	paddd	%xmm1, %xmm2
	pandn	%xmm0, %xmm1
	movdqa	%xmm2, %xmm3
	psrld	$31, %xmm2
	paddd	%xmm4, %xmm1
	pand	%xmm0, %xmm3
	paddd	%xmm3, %xmm2
	movups	%xmm2, 720(%r15)
	movdqu	976(%r15), %xmm2
	paddd	%xmm5, %xmm2
	movdqa	%xmm2, %xmm3
	psrld	$31, %xmm2
	pand	%xmm0, %xmm3
	paddd	%xmm3, %xmm2
	movups	%xmm2, 848(%r15)
	movdqa	%xmm1, %xmm2
	psrld	$31, %xmm1
	pand	%xmm0, %xmm2
	paddd	%xmm2, %xmm1
	movdqu	608(%r15), %xmm2
	movups	%xmm1, 976(%r15)
	movdqu	736(%r15), %xmm1
	movdqu	864(%r15), %xmm5
	movdqu	992(%r15), %xmm3
	movl	$5, 36(%rsp)
	movdqu	992(%r15), %xmm7
	movl	$4, 32(%rsp)
	pandn	%xmm0, %xmm3
	movl	$16, -80(%rsp)
	paddd	%xmm2, %xmm3
	paddd	%xmm7, %xmm2
	movdqu	1008(%r15), %xmm7
	movdqa	%xmm3, %xmm4
	pand	%xmm0, %xmm3
	psrld	$31, %xmm4
	paddd	%xmm4, %xmm3
	movups	%xmm3, 608(%r15)
	movdqa	%xmm5, %xmm3
	paddd	%xmm1, %xmm3
	pandn	%xmm0, %xmm1
	movdqa	%xmm3, %xmm4
	pand	%xmm0, %xmm3
	paddd	%xmm5, %xmm1
	movdqu	880(%r15), %xmm5
	psrld	$31, %xmm4
	paddd	%xmm4, %xmm3
	movups	%xmm3, 736(%r15)
	movdqa	%xmm2, %xmm3
	pand	%xmm0, %xmm2
	psrld	$31, %xmm3
	paddd	%xmm3, %xmm2
	movdqu	1008(%r15), %xmm3
	movups	%xmm2, 864(%r15)
	movdqa	%xmm1, %xmm2
	pand	%xmm0, %xmm1
	psrld	$31, %xmm2
	pandn	%xmm0, %xmm3
	paddd	%xmm2, %xmm1
	movdqu	624(%r15), %xmm2
	movups	%xmm1, 992(%r15)
	movdqu	752(%r15), %xmm1
	paddd	%xmm2, %xmm3
	paddd	%xmm7, %xmm2
	movdqa	%xmm3, %xmm4
	pand	%xmm0, %xmm3
	psrld	$31, %xmm4
	paddd	%xmm4, %xmm3
	movups	%xmm3, 624(%r15)
	movdqa	%xmm1, %xmm3
	pandn	%xmm0, %xmm1
	paddd	%xmm5, %xmm3
	paddd	%xmm5, %xmm1
	movdqa	%xmm3, %xmm4
	pand	%xmm0, %xmm3
	psrld	$31, %xmm4
	paddd	%xmm4, %xmm3
	movups	%xmm3, 752(%r15)
	movdqa	%xmm2, %xmm3
	pand	%xmm0, %xmm2
	psrld	$31, %xmm3
	paddd	%xmm3, %xmm2
	movups	%xmm2, 880(%r15)
	movdqa	%xmm1, %xmm2
	pand	%xmm0, %xmm1
	psrld	$31, %xmm2
	paddd	%xmm2, %xmm1
	movups	%xmm1, 1008(%r15)
	.p2align 4,,10
	.p2align 3
.L38:
	movl	-80(%rsp), %eax
	movl	$1, %ebp
	leal	(%rax,%rax), %ecx
	movl	%ecx, -56(%rsp)
	addl	%eax, %ecx
	movslq	-56(%rsp), %r15
	movl	%ecx, -76(%rsp)
	leal	0(,%rax,4), %ecx
	movl	%ecx, -52(%rsp)
	cmpl	$1, %eax
	je	.L3
	movslq	-80(%rsp), %rbp
	leaq	0(,%r15,4), %r8
	leaq	16(%r8), %r14
	leaq	(%r12,%r8), %rax
	leaq	0(,%rbp,4), %rdi
	leaq	(%r15,%rbp), %r9
	leaq	16(%rdi), %rbx
	cmpq	%rdi, %r14
	leaq	(%r12,%rdi), %rdx
	setle	%r11b
	cmpq	%rbx, %r8
	leaq	16(,%r9,4), %r13
	setge	%sil
	leaq	0(,%r9,4), %r10
	orl	%r11d, %esi
	movq	24(%rsp), %r11
	leaq	(%r12,%r10), %rcx
	cmpq	%r11, %rax
	setnb	%r11b
	andl	%esi, %r11d
	movq	24(%rsp), %rsi
	cmpq	%rsi, %rdx
	setnb	%sil
	andl	%r11d, %esi
	cmpq	%r8, %r13
	setle	%r8b
	cmpq	%r10, %r14
	setle	%r11b
	orl	%r11d, %r8d
	andl	%r8d, %esi
	cmpq	%rdi, %r13
	setle	%dil
	cmpq	%rbx, %r10
	setge	%r8b
	orl	%r8d, %edi
	testb	%dil, %sil
	je	.L3
	movq	24(%rsp), %rsi
	cmpq	%rsi, %rcx
	jb	.L3
	movl	-80(%rsp), %ebx
	leal	-1(%rbx), %edi
	movl	%ebx, %esi
	cmpl	$2, %edi
	jbe	.L49
	movl	%ebx, %edi
	xorl	%esi, %esi
	shrl	$2, %edi
	salq	$4, %rdi
	.p2align 4,,10
	.p2align 3
.L5:
	movdqu	(%r12,%rsi), %xmm6
	movdqu	(%rdx,%rsi), %xmm2
	movdqu	(%rax,%rsi), %xmm5
	movdqu	(%rcx,%rsi), %xmm1
	movdqa	%xmm6, %xmm3
	paddd	%xmm2, %xmm3
	pandn	%xmm0, %xmm2
	movdqa	%xmm3, %xmm4
	pand	%xmm0, %xmm3
	paddd	%xmm6, %xmm2
	psrld	$31, %xmm4
	paddd	%xmm4, %xmm3
	movups	%xmm3, (%r12,%rsi)
	movdqa	%xmm5, %xmm3
	paddd	%xmm1, %xmm3
	pandn	%xmm0, %xmm1
	movdqa	%xmm3, %xmm4
	pand	%xmm0, %xmm3
	paddd	%xmm5, %xmm1
	psrld	$31, %xmm4
	paddd	%xmm4, %xmm3
	movups	%xmm3, (%rdx,%rsi)
	movdqa	%xmm2, %xmm3
	pand	%xmm0, %xmm2
	psrld	$31, %xmm3
	paddd	%xmm3, %xmm2
	movups	%xmm2, (%rax,%rsi)
	movdqa	%xmm1, %xmm2
	pand	%xmm0, %xmm1
	psrld	$31, %xmm2
	paddd	%xmm2, %xmm1
	movups	%xmm1, (%rcx,%rsi)
	addq	$16, %rsi
	cmpq	%rsi, %rdi
	jne	.L5
	movl	-80(%rsp), %ecx
	movl	%ecx, %eax
	andl	$-4, %eax
	movl	%eax, %edx
	cmpl	%eax, %ecx
	je	.L6
	subl	%eax, %ecx
	movl	%ecx, %esi
	cmpl	$1, %ecx
	je	.L7
.L4:
	leaq	(%r12,%rdx,4), %r8
	leaq	0(%rbp,%rdx), %rcx
	addq	%rdx, %r9
	movq	.LC1(%rip), %xmm3
	movq	(%r8), %xmm7
	leaq	(%r12,%rcx,4), %rdi
	leaq	(%r15,%rdx), %rcx
	movq	(%rdi), %xmm2
	leaq	(%r12,%rcx,4), %rcx
	leaq	(%r12,%r9,4), %rdx
	movdqa	%xmm7, %xmm4
	movq	(%rcx), %xmm5
	movq	(%rdx), %xmm1
	paddd	%xmm2, %xmm4
	pandn	%xmm3, %xmm2
	paddd	%xmm7, %xmm2
	movdqa	%xmm4, %xmm6
	pand	%xmm3, %xmm4
	psrld	$31, %xmm6
	paddd	%xmm6, %xmm4
	movq	%xmm4, (%r8)
	movdqa	%xmm5, %xmm4
	paddd	%xmm1, %xmm4
	pandn	%xmm3, %xmm1
	paddd	%xmm5, %xmm1
	movdqa	%xmm4, %xmm6
	pand	%xmm3, %xmm4
	psrld	$31, %xmm6
	paddd	%xmm6, %xmm4
	movq	%xmm4, (%rdi)
	movdqa	%xmm2, %xmm4
	pand	%xmm3, %xmm2
	psrld	$31, %xmm4
	paddd	%xmm4, %xmm2
	movq	%xmm2, (%rcx)
	movdqa	%xmm1, %xmm2
	pand	%xmm3, %xmm1
	psrld	$31, %xmm2
	paddd	%xmm2, %xmm1
	movq	%xmm1, (%rdx)
	testb	$1, %sil
	je	.L6
	andl	$-2, %esi
	addl	%esi, %eax
.L7:
	movl	-80(%rsp), %ecx
	movslq	%eax, %rdx
	leaq	(%r12,%rdx,4), %rbx
	leal	(%rcx,%rax), %edx
	movl	-56(%rsp), %ecx
	movl	(%rbx), %r10d
	movslq	%edx, %rdx
	addl	%eax, %ecx
	leaq	(%r12,%rdx,4), %r11
	movslq	%ecx, %rcx
	movl	(%r11), %edx
	leaq	(%r12,%rcx,4), %r9
	movl	-76(%rsp), %ecx
	movl	(%r9), %r8d
	addl	%ecx, %eax
	leal	(%r10,%rdx), %ecx
	notl	%edx
	cltq
	movl	%ecx, %esi
	andl	$2147483647, %ecx
	andl	$2147483647, %edx
	leaq	(%r12,%rax,4), %rdi
	shrl	$31, %esi
	addl	%r10d, %edx
	movl	(%rdi), %eax
	addl	%ecx, %esi
	movl	%esi, (%rbx)
	leal	(%r8,%rax), %ecx
	notl	%eax
	movl	%ecx, %esi
	andl	$2147483647, %ecx
	andl	$2147483647, %eax
	shrl	$31, %esi
	addl	%r8d, %eax
	addl	%ecx, %esi
	movl	%edx, %ecx
	andl	$2147483647, %edx
	shrl	$31, %ecx
	movl	%esi, (%r11)
	addl	%edx, %ecx
	movl	%eax, %edx
	andl	$2147483647, %eax
	shrl	$31, %edx
	movl	%ecx, (%r9)
	addl	%edx, %eax
	movl	%eax, (%rdi)
.L6:
	movslq	-52(%rsp), %r14
	leaq	(%r15,%r14), %r9
	leaq	0(%rbp,%r14), %rcx
	leaq	(%r9,%rbp), %r13
	leaq	0(,%rcx,4), %rdx
	movq	%rcx, -104(%rsp)
	leaq	16(,%r13,4), %rax
	leaq	16(%rdx), %r10
	movq	%rax, -120(%rsp)
	leaq	0(,%r9,4), %rax
	leaq	0(,%r14,4), %rcx
	leaq	16(%rax), %r11
	leaq	16(%rcx), %rbx
	cmpq	%rdx, %r11
	leaq	0(,%r13,4), %rsi
	setle	%dil
	cmpq	%r10, %rax
	setge	%r8b
	orl	%r8d, %edi
	cmpq	%rcx, %r11
	setle	%r8b
	cmpq	%rbx, %rax
	setge	-112(%rsp)
	orb	-112(%rsp), %r8b
	andl	%edi, %r8d
	cmpq	%rcx, %r10
	setle	%dil
	cmpq	%rbx, %rdx
	setge	-112(%rsp)
	orb	-112(%rsp), %dil
	andl	%r8d, %edi
	cmpq	%rax, -120(%rsp)
	setle	%r8b
	cmpq	%rsi, %r11
	setle	%r11b
	orl	%r11d, %r8d
	andl	%edi, %r8d
	cmpq	%rdx, -120(%rsp)
	setle	%dil
	cmpq	%r10, %rsi
	setge	%r10b
	orl	%r10d, %edi
	testb	%dil, %r8b
	je	.L137
	cmpq	%rcx, -120(%rsp)
	setle	%dil
	cmpq	%rbx, %rsi
	setge	%r8b
	orb	%dil, %r8b
	je	.L137
	movl	-80(%rsp), %ebx
	leal	-1(%rbx), %r8d
	movl	%ebx, %edi
	cmpl	$2, %r8d
	jbe	.L50
	shrl	$2, %ebx
	addq	%r12, %rcx
	addq	%r12, %rdx
	addq	%r12, %rax
	movl	%ebx, %r8d
	addq	%r12, %rsi
	xorl	%edi, %edi
	salq	$4, %r8
	.p2align 4,,10
	.p2align 3
.L13:
	movdqu	(%rsi,%rdi), %xmm6
	movdqu	(%rcx,%rdi), %xmm2
	movdqu	(%rdx,%rdi), %xmm1
	movdqu	(%rax,%rdi), %xmm5
	movdqa	%xmm6, %xmm3
	pandn	%xmm0, %xmm3
	paddd	%xmm2, %xmm3
	paddd	%xmm6, %xmm2
	movdqa	%xmm3, %xmm4
	pand	%xmm0, %xmm3
	psrld	$31, %xmm4
	paddd	%xmm4, %xmm3
	movups	%xmm3, (%rcx,%rdi)
	movdqa	%xmm1, %xmm3
	pandn	%xmm0, %xmm1
	paddd	%xmm5, %xmm3
	paddd	%xmm5, %xmm1
	movdqa	%xmm3, %xmm4
	pand	%xmm0, %xmm3
	psrld	$31, %xmm4
	paddd	%xmm4, %xmm3
	movups	%xmm3, (%rdx,%rdi)
	movdqa	%xmm2, %xmm3
	pand	%xmm0, %xmm2
	psrld	$31, %xmm3
	paddd	%xmm3, %xmm2
	movups	%xmm2, (%rax,%rdi)
	movdqa	%xmm1, %xmm2
	pand	%xmm0, %xmm1
	psrld	$31, %xmm2
	paddd	%xmm2, %xmm1
	movups	%xmm1, (%rsi,%rdi)
	addq	$16, %rdi
	cmpq	%r8, %rdi
	jne	.L13
	movl	-80(%rsp), %ecx
	movl	%ecx, %edx
	andl	$-4, %edx
	movl	%edx, %eax
	cmpl	%edx, %ecx
	je	.L138
	movl	%ecx, %edi
	subl	%edx, %edi
	cmpl	$1, %edi
	je	.L16
.L12:
	leaq	(%r14,%rax), %rcx
	addq	%rax, %r13
	addq	%rax, %r9
	movq	.LC1(%rip), %xmm3
	leaq	(%r12,%rcx,4), %rsi
	movq	-104(%rsp), %rcx
	leaq	(%r12,%r9,4), %r8
	movq	(%rsi), %xmm2
	movq	(%r8), %xmm5
	addq	%rax, %rcx
	leaq	(%r12,%r13,4), %rax
	movq	(%rax), %xmm7
	leaq	(%r12,%rcx,4), %rcx
	movq	(%rcx), %xmm1
	movdqa	%xmm7, %xmm4
	pandn	%xmm3, %xmm4
	paddd	%xmm2, %xmm4
	paddd	%xmm7, %xmm2
	movdqa	%xmm4, %xmm6
	pand	%xmm3, %xmm4
	psrld	$31, %xmm6
	paddd	%xmm6, %xmm4
	movq	%xmm4, (%rsi)
	movdqa	%xmm1, %xmm4
	pandn	%xmm3, %xmm1
	paddd	%xmm5, %xmm4
	paddd	%xmm5, %xmm1
	movdqa	%xmm4, %xmm6
	pand	%xmm3, %xmm4
	psrld	$31, %xmm6
	paddd	%xmm6, %xmm4
	movq	%xmm4, (%rcx)
	movdqa	%xmm2, %xmm4
	pand	%xmm3, %xmm2
	psrld	$31, %xmm4
	paddd	%xmm4, %xmm2
	movq	%xmm2, (%r8)
	movdqa	%xmm1, %xmm2
	pand	%xmm3, %xmm1
	psrld	$31, %xmm2
	paddd	%xmm2, %xmm1
	movq	%xmm1, (%rax)
	testb	$1, %dil
	je	.L138
	andl	$-2, %edi
	addl	%edi, %edx
.L16:
	movl	-52(%rsp), %eax
	movl	-56(%rsp), %edi
	addl	%eax, %edx
	movslq	%edx, %rax
	leal	(%rdi,%rdx), %esi
	movl	-76(%rsp), %edi
	leaq	(%r12,%rax,4), %rbx
	movl	-80(%rsp), %eax
	movslq	%esi, %rsi
	movl	(%rbx), %ecx
	leaq	(%r12,%rsi,4), %r9
	addl	%edx, %eax
	addl	%edi, %edx
	movl	(%r9), %r8d
	movslq	%edx, %rdx
	cltq
	leaq	(%r12,%rdx,4), %rdi
	leaq	(%r12,%rax,4), %r11
	movl	(%rdi), %r10d
	movl	(%r11), %eax
	movl	%r10d, %edx
	notl	%edx
	andl	$2147483647, %edx
	addl	%ecx, %edx
	addl	%r10d, %ecx
	movl	%edx, %esi
	andl	$2147483647, %edx
	shrl	$31, %esi
	addl	%edx, %esi
	leal	(%rax,%r8), %edx
	notl	%eax
	movl	%esi, (%rbx)
	movl	%edx, %esi
	andl	$2147483647, %edx
	andl	$2147483647, %eax
	shrl	$31, %esi
	addl	%r8d, %eax
	addl	%edx, %esi
	movl	%ecx, %edx
	andl	$2147483647, %ecx
	shrl	$31, %edx
	movl	%esi, (%r11)
	addl	%ecx, %edx
	movl	%edx, (%r9)
	movl	%eax, %edx
	andl	$2147483647, %eax
	shrl	$31, %edx
	addl	%edx, %eax
	movl	%eax, (%rdi)
.L138:
	movl	-52(%rsp), %eax
	leaq	(%r14,%r14), %r8
	addl	%eax, %eax
	movl	%eax, -120(%rsp)
	leaq	(%r15,%rbp), %rax
	leaq	(%rax,%r8), %rcx
	leaq	(%r15,%r8), %rax
	leaq	16(,%rcx,4), %r10
	movq	%rax, -96(%rsp)
	salq	$2, %rax
	leaq	0(,%rcx,4), %rsi
	leaq	16(%rax), %r9
	cmpq	%rax, %r10
	movq	%rcx, -104(%rsp)
	leaq	0(%rbp,%r8), %rcx
	setle	%r11b
	cmpq	%rsi, %r9
	movq	%rcx, -88(%rsp)
	leaq	0(,%rcx,4), %rdx
	setle	%dil
	leaq	16(%rdx), %r13
	leaq	0(,%r14,8), %rcx
	orl	%edi, %r11d
	cmpq	%rdx, %r10
	leaq	16(%rcx), %rbx
	setle	%dil
	cmpq	%r13, %rsi
	setge	-112(%rsp)
	orb	-112(%rsp), %dil
	andl	%r11d, %edi
	cmpq	%rcx, %r10
	setle	%r10b
	cmpq	%rbx, %rsi
	setge	%r11b
	orl	%r11d, %r10d
	andl	%edi, %r10d
	cmpq	%rdx, %r9
	setle	%dil
	cmpq	%r13, %rax
	setge	%r11b
	orl	%r11d, %edi
	andl	%r10d, %edi
	cmpq	%rcx, %r9
	setle	%r9b
	cmpq	%rbx, %rax
	setge	%r10b
	orl	%r10d, %r9d
	testb	%r9b, %dil
	je	.L139
	cmpq	%rcx, %r13
	setle	%dil
	cmpq	%rbx, %rdx
	setge	%r9b
	orb	%dil, %r9b
	je	.L139
	movl	-80(%rsp), %ebx
	leal	-1(%rbx), %edi
	movl	%ebx, %r9d
	cmpl	$2, %edi
	jbe	.L52
	shrl	$2, %ebx
	addq	%r12, %rcx
	addq	%r12, %rdx
	addq	%r12, %rax
	movl	%ebx, %r9d
	addq	%r12, %rsi
	xorl	%edi, %edi
	salq	$4, %r9
	.p2align 4,,10
	.p2align 3
.L22:
	movdqu	(%rdx,%rdi), %xmm3
	movdqu	(%rsi,%rdi), %xmm2
	movdqu	(%rcx,%rdi), %xmm6
	movdqu	(%rax,%rdi), %xmm4
	movdqa	%xmm3, %xmm1
	psrld	$16, %xmm3
	pslld	$15, %xmm1
	pand	.LC3(%rip), %xmm1
	por	%xmm3, %xmm1
	movdqa	%xmm2, %xmm3
	pslld	$15, %xmm3
	pand	.LC3(%rip), %xmm3
	psrld	$16, %xmm2
	por	%xmm2, %xmm3
	movdqa	%xmm3, %xmm5
	pxor	%xmm0, %xmm5
	paddd	%xmm1, %xmm5
	paddd	%xmm3, %xmm1
	movdqa	%xmm5, %xmm2
	movdqa	%xmm1, %xmm3
	pand	%xmm0, %xmm5
	psrld	$31, %xmm3
	psrld	$31, %xmm2
	pand	%xmm0, %xmm1
	paddd	%xmm5, %xmm2
	paddd	%xmm3, %xmm1
	movdqa	%xmm6, %xmm3
	paddd	%xmm2, %xmm3
	pandn	%xmm0, %xmm2
	movdqa	%xmm3, %xmm5
	pand	%xmm0, %xmm3
	paddd	%xmm6, %xmm2
	psrld	$31, %xmm5
	paddd	%xmm5, %xmm3
	movups	%xmm3, (%rcx,%rdi)
	movdqa	%xmm4, %xmm3
	paddd	%xmm1, %xmm3
	pandn	%xmm0, %xmm1
	movdqa	%xmm3, %xmm5
	pand	%xmm0, %xmm3
	paddd	%xmm4, %xmm1
	psrld	$31, %xmm5
	paddd	%xmm5, %xmm3
	movups	%xmm3, (%rdx,%rdi)
	movdqa	%xmm2, %xmm3
	pand	%xmm0, %xmm2
	psrld	$31, %xmm3
	paddd	%xmm3, %xmm2
	movups	%xmm2, (%rax,%rdi)
	movdqa	%xmm1, %xmm2
	pand	%xmm0, %xmm1
	psrld	$31, %xmm2
	paddd	%xmm2, %xmm1
	movups	%xmm1, (%rsi,%rdi)
	addq	$16, %rdi
	cmpq	%r9, %rdi
	jne	.L22
	movl	-80(%rsp), %ecx
	movl	%ecx, %eax
	andl	$-4, %eax
	movl	%eax, %edx
	cmpl	%eax, %ecx
	je	.L26
	movl	%ecx, %r9d
	subl	%eax, %r9d
	cmpl	$1, %r9d
	je	.L25
.L21:
	movq	-88(%rsp), %rcx
	movq	-104(%rsp), %r13
	addq	%rdx, %r8
	movq	.LC4(%rip), %xmm7
	leaq	(%r12,%r8,4), %rdi
	movq	.LC1(%rip), %xmm5
	addq	%rdx, %rcx
	addq	%rdx, %r13
	movq	(%rdi), %xmm6
	leaq	(%r12,%rcx,4), %rsi
	movq	-96(%rsp), %rcx
	movq	(%rsi), %xmm3
	addq	%rdx, %rcx
	leaq	(%r12,%r13,4), %rdx
	movdqa	%xmm3, %xmm1
	movq	(%rdx), %xmm2
	psrld	$16, %xmm3
	leaq	(%r12,%rcx,4), %rcx
	pslld	$15, %xmm1
	movq	(%rcx), %xmm4
	pand	%xmm7, %xmm1
	por	%xmm3, %xmm1
	movdqa	%xmm2, %xmm3
	pslld	$15, %xmm3
	psrld	$16, %xmm2
	pand	%xmm7, %xmm3
	movq	.LC1(%rip), %xmm7
	por	%xmm2, %xmm3
	pxor	%xmm3, %xmm5
	paddd	%xmm1, %xmm5
	paddd	%xmm3, %xmm1
	movdqa	%xmm5, %xmm2
	pand	%xmm7, %xmm5
	movdqa	%xmm1, %xmm3
	psrld	$31, %xmm2
	psrld	$31, %xmm3
	pand	%xmm7, %xmm1
	paddd	%xmm5, %xmm2
	paddd	%xmm3, %xmm1
	movdqa	%xmm6, %xmm3
	paddd	%xmm2, %xmm3
	pandn	%xmm7, %xmm2
	paddd	%xmm6, %xmm2
	movdqa	%xmm3, %xmm5
	pand	%xmm7, %xmm3
	psrld	$31, %xmm5
	paddd	%xmm5, %xmm3
	movq	%xmm3, (%rdi)
	movdqa	%xmm4, %xmm3
	paddd	%xmm1, %xmm3
	pandn	%xmm7, %xmm1
	paddd	%xmm4, %xmm1
	movdqa	%xmm3, %xmm5
	pand	%xmm7, %xmm3
	psrld	$31, %xmm5
	paddd	%xmm5, %xmm3
	movq	%xmm3, (%rsi)
	movdqa	%xmm2, %xmm3
	pand	%xmm7, %xmm2
	psrld	$31, %xmm3
	paddd	%xmm3, %xmm2
	movq	%xmm2, (%rcx)
	movdqa	%xmm1, %xmm2
	pand	%xmm7, %xmm1
	psrld	$31, %xmm2
	paddd	%xmm2, %xmm1
	movq	%xmm1, (%rdx)
	testb	$1, %r9b
	je	.L26
	andl	$-2, %r9d
	addl	%r9d, %eax
.L25:
	movl	-120(%rsp), %ecx
	addl	%ecx, %eax
	movl	-80(%rsp), %ecx
	movslq	%eax, %rdx
	leaq	(%r12,%rdx,4), %r10
	leal	(%rcx,%rax), %edx
	movl	-56(%rsp), %ecx
	movslq	%edx, %rdx
	movl	(%r10), %r8d
	leaq	(%r12,%rdx,4), %r9
	leal	(%rcx,%rax), %edx
	movl	-76(%rsp), %ecx
	movl	(%r9), %r11d
	movslq	%edx, %rdx
	addl	%ecx, %eax
	leaq	(%r12,%rdx,4), %rdi
	cltq
	movl	%r11d, %edx
	shrl	$16, %r11d
	movl	(%rdi), %esi
	leaq	(%r12,%rax,4), %rcx
	sall	$15, %edx
	movl	(%rcx), %eax
	andl	$2147450880, %edx
	orl	%r11d, %edx
	movl	%eax, %r11d
	shrl	$16, %eax
	sall	$15, %r11d
	andl	$2147450880, %r11d
	orl	%eax, %r11d
	movl	%r11d, %eax
	xorl	$2147483647, %eax
	addl	%edx, %eax
	addl	%r11d, %edx
	movl	%eax, %ebx
	movl	%edx, %r11d
	andl	$2147483647, %eax
	andl	$2147483647, %edx
	shrl	$31, %ebx
	shrl	$31, %r11d
	addl	%ebx, %eax
	addl	%r11d, %edx
	leal	(%r8,%rax), %r11d
	notl	%eax
	movl	%r11d, %ebx
	andl	$2147483647, %eax
	andl	$2147483647, %r11d
	shrl	$31, %ebx
	addl	%r8d, %eax
	addl	%r11d, %ebx
	movl	%eax, %r8d
	andl	$2147483647, %eax
	movl	%ebx, (%r10)
	leal	(%rsi,%rdx), %r10d
	notl	%edx
	shrl	$31, %r8d
	andl	$2147483647, %edx
	addl	%eax, %r8d
	movl	%r10d, %r11d
	andl	$2147483647, %r10d
	leal	(%rdx,%rsi), %eax
	shrl	$31, %r11d
	movl	%eax, %edx
	addl	%r10d, %r11d
	shrl	$31, %edx
	andl	$2147483647, %eax
	movl	%r11d, (%r9)
	addl	%eax, %edx
	movl	%r8d, (%rdi)
	movl	%edx, (%rcx)
.L26:
	movl	-120(%rsp), %eax
	movl	-52(%rsp), %ecx
	addl	%ecx, %eax
	movl	%eax, -96(%rsp)
.L24:
	leaq	(%r14,%r14,2), %rbx
	leaq	0(%rbp,%r15), %rax
	leaq	(%rbx,%rbp), %rcx
	leaq	(%rax,%rbx), %r13
	leaq	0(,%rcx,4), %rdx
	leaq	0(,%r13,4), %rsi
	movq	%rcx, -104(%rsp)
	leaq	(%rbx,%r15), %rax
	leaq	16(%rdx), %r11
	leaq	16(,%r13,4), %r9
	leaq	0(,%rbx,4), %rcx
	movq	%rax, -88(%rsp)
	salq	$2, %rax
	cmpq	%rsi, %r11
	leaq	16(%rcx), %rdi
	leaq	16(%rax), %r8
	setle	%r10b
	cmpq	%r9, %rdx
	movq	%rdi, -120(%rsp)
	setge	%dil
	orl	%edi, %r10d
	cmpq	%r8, %rsi
	setge	%dil
	cmpq	%r9, %rax
	setge	-112(%rsp)
	orb	-112(%rsp), %dil
	andl	%r10d, %edi
	cmpq	%rsi, -120(%rsp)
	setle	%r10b
	cmpq	%r9, %rcx
	setge	%r9b
	orl	%r10d, %r9d
	andl	%edi, %r9d
	cmpq	%rax, %r11
	setle	%dil
	cmpq	%r8, %rdx
	setge	%r10b
	orl	%r10d, %edi
	andl	%r9d, %edi
	cmpq	%rax, -120(%rsp)
	setle	%r9b
	cmpq	%r8, %rcx
	setge	%r8b
	orl	%r8d, %r9d
	testb	%r9b, %dil
	je	.L140
	cmpq	%rdx, -120(%rsp)
	setle	%dil
	cmpq	%r11, %rcx
	setge	%r8b
	orb	%dil, %r8b
	je	.L140
	movl	-80(%rsp), %r9d
	leal	-1(%r9), %edi
	movl	%r9d, %r8d
	cmpl	$2, %edi
	jbe	.L53
	shrl	$2, %r9d
	addq	%r12, %rcx
	addq	%r12, %rdx
	addq	%r12, %rax
	movl	%r9d, %r8d
	addq	%r12, %rsi
	xorl	%edi, %edi
	salq	$4, %r8
	.p2align 4,,10
	.p2align 3
.L31:
	movdqu	(%rdx,%rdi), %xmm1
	movdqu	(%rsi,%rdi), %xmm2
	movdqu	(%rcx,%rdi), %xmm5
	movdqu	(%rax,%rdi), %xmm4
	movdqa	%xmm1, %xmm3
	psrld	$16, %xmm1
	pslld	$15, %xmm3
	pand	.LC3(%rip), %xmm3
	por	%xmm1, %xmm3
	movdqa	%xmm2, %xmm1
	pslld	$15, %xmm1
	pand	.LC3(%rip), %xmm1
	psrld	$16, %xmm2
	por	%xmm2, %xmm1
	movdqa	%xmm1, %xmm6
	pxor	%xmm0, %xmm1
	paddd	%xmm3, %xmm6
	paddd	%xmm3, %xmm1
	movdqa	%xmm6, %xmm2
	psrld	$31, %xmm6
	movdqa	%xmm1, %xmm3
	pand	%xmm0, %xmm2
	psrld	$31, %xmm1
	pand	%xmm0, %xmm3
	paddd	%xmm6, %xmm2
	paddd	%xmm3, %xmm1
	movdqa	%xmm2, %xmm3
	paddd	%xmm5, %xmm2
	pandn	%xmm0, %xmm3
	paddd	%xmm5, %xmm3
	movdqa	%xmm3, %xmm6
	psrld	$31, %xmm3
	pand	%xmm0, %xmm6
	paddd	%xmm6, %xmm3
	movups	%xmm3, (%rcx,%rdi)
	movdqa	%xmm1, %xmm3
	pandn	%xmm0, %xmm1
	paddd	%xmm4, %xmm3
	paddd	%xmm4, %xmm1
	movdqa	%xmm3, %xmm6
	psrld	$31, %xmm3
	pand	%xmm0, %xmm6
	paddd	%xmm6, %xmm3
	movups	%xmm3, (%rdx,%rdi)
	movdqa	%xmm2, %xmm3
	psrld	$31, %xmm2
	pand	%xmm0, %xmm3
	paddd	%xmm3, %xmm2
	movups	%xmm2, (%rax,%rdi)
	movdqa	%xmm1, %xmm2
	psrld	$31, %xmm1
	pand	%xmm0, %xmm2
	paddd	%xmm2, %xmm1
	movups	%xmm1, (%rsi,%rdi)
	addq	$16, %rdi
	cmpq	%r8, %rdi
	jne	.L31
	movl	-80(%rsp), %ecx
	movl	%ecx, %edx
	andl	$-4, %edx
	movl	%edx, %eax
	cmpl	%ecx, %edx
	je	.L36
	movl	%ecx, %r8d
	subl	%edx, %r8d
	cmpl	$1, %r8d
	je	.L33
.L30:
	movq	.LC4(%rip), %xmm7
	leaq	(%rax,%rbx), %rcx
	leaq	(%r12,%rcx,4), %rdi
	movq	-104(%rsp), %rcx
	movq	(%rdi), %xmm5
	addq	%rax, %rcx
	leaq	(%r12,%rcx,4), %rsi
	movq	-88(%rsp), %rcx
	movq	(%rsi), %xmm1
	addq	%rax, %rcx
	addq	%r13, %rax
	leaq	(%r12,%rax,4), %rax
	movdqa	%xmm1, %xmm3
	leaq	(%r12,%rcx,4), %rcx
	movq	(%rax), %xmm2
	pslld	$15, %xmm3
	psrld	$16, %xmm1
	movq	(%rcx), %xmm4
	pand	%xmm7, %xmm3
	por	%xmm1, %xmm3
	movdqa	%xmm2, %xmm1
	pslld	$15, %xmm1
	psrld	$16, %xmm2
	pand	%xmm7, %xmm1
	movq	.LC1(%rip), %xmm7
	por	%xmm2, %xmm1
	movq	.LC1(%rip), %xmm2
	movdqa	%xmm1, %xmm6
	pxor	%xmm7, %xmm1
	paddd	%xmm3, %xmm6
	paddd	%xmm3, %xmm1
	movdqa	%xmm7, %xmm3
	pand	%xmm6, %xmm2
	psrld	$31, %xmm6
	pand	%xmm1, %xmm3
	paddd	%xmm6, %xmm2
	psrld	$31, %xmm1
	movdqa	%xmm7, %xmm6
	paddd	%xmm3, %xmm1
	movdqa	%xmm2, %xmm3
	paddd	%xmm5, %xmm2
	pandn	%xmm7, %xmm3
	paddd	%xmm5, %xmm3
	pand	%xmm3, %xmm6
	psrld	$31, %xmm3
	paddd	%xmm6, %xmm3
	movdqa	%xmm7, %xmm6
	movq	%xmm3, (%rdi)
	movdqa	%xmm1, %xmm3
	pandn	%xmm7, %xmm1
	paddd	%xmm4, %xmm3
	paddd	%xmm4, %xmm1
	pand	%xmm3, %xmm6
	psrld	$31, %xmm3
	paddd	%xmm6, %xmm3
	movq	%xmm3, (%rsi)
	movdqa	%xmm7, %xmm3
	pand	%xmm2, %xmm3
	psrld	$31, %xmm2
	paddd	%xmm3, %xmm2
	movq	%xmm2, (%rcx)
	movdqa	%xmm1, %xmm2
	pand	%xmm7, %xmm1
	psrld	$31, %xmm2
	paddd	%xmm2, %xmm1
	movq	%xmm1, (%rax)
	testb	$1, %r8b
	je	.L36
	andl	$-2, %r8d
	addl	%r8d, %edx
.L33:
	movl	-96(%rsp), %eax
	movl	-56(%rsp), %ecx
	addl	%eax, %edx
	movslq	%edx, %rax
	addl	%edx, %ecx
	leaq	(%r12,%rax,4), %r11
	movl	-80(%rsp), %eax
	movslq	%ecx, %rcx
	leaq	(%r12,%rcx,4), %r8
	movl	-76(%rsp), %ecx
	movl	(%r11), %r9d
	addl	%edx, %eax
	movl	(%r8), %edi
	cltq
	addl	%ecx, %edx
	leaq	(%r12,%rax,4), %r10
	movslq	%edx, %rdx
	movl	(%r10), %eax
	leaq	(%r12,%rdx,4), %rsi
	movl	(%rsi), %edx
	movl	%eax, %ecx
	shrl	$16, %eax
	sall	$15, %ecx
	andl	$2147450880, %ecx
	orl	%eax, %ecx
	movl	%edx, %eax
	shrl	$16, %edx
	sall	$15, %eax
	andl	$2147450880, %eax
	orl	%edx, %eax
	leal	(%rax,%rcx), %ebx
	xorl	$2147483647, %eax
	movl	%ebx, %edx
	addl	%eax, %ecx
	shrl	$31, %ebx
	andl	$2147483647, %edx
	addl	%ebx, %edx
	movl	%ecx, %ebx
	shrl	$31, %ecx
	andl	$2147483647, %ebx
	leal	(%rbx,%rcx), %eax
	movl	%edx, %ecx
	addl	%r9d, %edx
	notl	%ecx
	andl	$2147483647, %ecx
	addl	%r9d, %ecx
	movl	%ecx, %ebx
	shrl	$31, %ecx
	andl	$2147483647, %ebx
	addl	%ecx, %ebx
	leal	(%rax,%rdi), %ecx
	notl	%eax
	movl	%ebx, (%r11)
	movl	%ecx, %r11d
	shrl	$31, %ecx
	andl	$2147483647, %eax
	andl	$2147483647, %r11d
	addl	%ecx, %r11d
	movl	%edx, %ecx
	shrl	$31, %edx
	andl	$2147483647, %ecx
	movl	%r11d, (%r10)
	addl	%edx, %ecx
	addl	%edi, %eax
	movl	%eax, %edx
	shrl	$31, %eax
	movl	%ecx, (%r8)
	andl	$2147483647, %edx
	addl	%eax, %edx
	movl	%edx, (%rsi)
.L36:
	cmpl	$4, 32(%rsp)
	jle	.L48
	movslq	-96(%rsp), %rdx
	leaq	4(%r14,%rbp), %rcx
	movq	%r12, -64(%rsp)
	leaq	32+omegas128(%rip), %rax
	movq	%rax, -104(%rsp)
	leaq	0(,%r14,4), %rax
	movslq	-76(%rsp), %r13
	movq	%rax, -24(%rsp)
	leaq	(%rdx,%r14), %rax
	addq	%rdx, %rcx
	leaq	(%r12,%rax,4), %rbx
	leaq	0(%rbp,%r14), %rax
	salq	$2, %rcx
	movq	%rbp, -88(%rsp)
	addq	%rdx, %rax
	leaq	4(%rdx,%r14), %rdx
	movq	%rcx, -112(%rsp)
	movq	%rbx, %r14
	leaq	0(,%rdx,4), %rcx
	salq	$2, %rax
	movq	%rcx, -120(%rsp)
	movl	32(%rsp), %ecx
	movq	%rax, -72(%rsp)
	addq	%r12, %rax
	leal	-5(%rcx), %edx
	leaq	40+omegas128(%rip), %rcx
	leaq	(%rcx,%rdx,8), %rcx
	movq	%rcx, -32(%rsp)
	movl	-80(%rsp), %ecx
	leal	-1(%rcx), %edi
	movl	%ecx, %edx
	movl	%edi, -8(%rsp)
	leaq	0(,%r15,4), %rdi
	shrl	$2, %edx
	movq	%rdi, -48(%rsp)
	subq	%r12, %rdi
	salq	$4, %rdx
	movq	%rdi, -40(%rsp)
	movl	%ecx, %edi
	andl	$3, %ecx
	andl	$-4, %edi
	movq	%rdx, (%rsp)
	movl	%edi, 16(%rsp)
	leal	1(%rdi), %edx
	addl	$2, %edi
	movl	%edx, 12(%rsp)
	movl	%edi, 20(%rsp)
	movl	%ecx, 8(%rsp)
	.p2align 4,,10
	.p2align 3
.L47:
	movq	-104(%rsp), %rcx
	movl	-52(%rsp), %edi
	addl	%edi, -96(%rsp)
	cmpl	$2, -8(%rsp)
	movl	(%rcx), %edx
	movl	4(%rcx), %ecx
	movq	%rdx, %xmm5
	movq	%rcx, %xmm4
	jbe	.L54
	movq	-40(%rsp), %rdi
	movq	-48(%rsp), %rbx
	movq	-112(%rsp), %r12
	movq	-120(%rsp), %r11
	leaq	(%rdi,%rax), %r10
	movq	-72(%rsp), %rbp
	leaq	(%rdi,%r14), %r8
	leaq	(%rbx,%r12), %r9
	cmpq	%r10, %r12
	leaq	(%rbx,%r11), %rdi
	leaq	-16(%r11), %rbx
	setle	%r11b
	cmpq	%r9, %rbp
	setge	%sil
	orl	%esi, %r11d
	cmpq	%r10, %rdi
	setle	%sil
	cmpq	%r8, %r9
	setle	%r12b
	orl	%r12d, %esi
	movq	-112(%rsp), %r12
	andl	%r11d, %esi
	movq	-120(%rsp), %r11
	cmpq	%r10, %r11
	setle	%r10b
	cmpq	%r9, %rbx
	setge	%r9b
	orl	%r9d, %r10d
	andl	%r10d, %esi
	cmpq	%r8, %r12
	setle	%r9b
	cmpq	%rdi, %rbp
	setge	%r10b
	orl	%r10d, %r9d
	andl	%esi, %r9d
	cmpq	%r8, %r11
	setle	%sil
	cmpq	%rdi, %rbx
	setge	%dil
	orl	%edi, %esi
	testb	%sil, %r9b
	je	.L54
	cmpq	%rbp, %r11
	setle	%sil
	cmpq	%rbx, %r12
	setle	%dil
	orb	%sil, %dil
	je	.L54
	pshufd	$0, %xmm5, %xmm5
	pshufd	$0, %xmm4, %xmm4
	movq	(%rsp), %r9
	xorl	%esi, %esi
	movq	-48(%rsp), %rdi
	movdqa	%xmm5, %xmm7
	movdqa	%xmm4, %xmm6
	movq	-88(%rsp), %rbp
	punpckldq	%xmm5, %xmm7
	punpckldq	%xmm4, %xmm6
	punpckhdq	%xmm5, %xmm5
	leaq	(%rdi,%r14), %r8
	punpckhdq	%xmm4, %xmm4
	addq	%rax, %rdi
	.p2align 4,,10
	.p2align 3
.L44:
	movdqu	(%rax,%rsi), %xmm8
	movdqa	%xmm7, %xmm1
	movdqa	%xmm5, %xmm10
	movdqu	(%rdi,%rsi), %xmm3
	movdqa	%xmm4, %xmm14
	movdqu	(%r14,%rsi), %xmm11
	movdqu	(%r8,%rsi), %xmm9
	movdqa	%xmm8, %xmm12
	punpckldq	%xmm8, %xmm12
	punpckhdq	%xmm8, %xmm8
	pmuludq	%xmm8, %xmm10
	pmuludq	%xmm12, %xmm1
	pmuludq	%xmm6, %xmm12
	movdqa	%xmm1, %xmm2
	movdqa	%xmm10, %xmm13
	shufps	$136, %xmm10, %xmm1
	movdqa	%xmm3, %xmm10
	psrlq	$31, %xmm13
	pand	%xmm0, %xmm1
	punpckldq	%xmm3, %xmm10
	psrlq	$31, %xmm2
	punpckhdq	%xmm3, %xmm3
	shufps	$136, %xmm13, %xmm2
	paddd	%xmm1, %xmm2
	movdqa	%xmm6, %xmm1
	pmuludq	%xmm10, %xmm1
	pmuludq	%xmm3, %xmm14
	pmuludq	%xmm7, %xmm10
	pmuludq	%xmm5, %xmm3
	movdqa	%xmm1, %xmm13
	movdqa	%xmm14, %xmm15
	shufps	$136, %xmm14, %xmm1
	pand	%xmm0, %xmm1
	psrlq	$31, %xmm13
	psrlq	$31, %xmm15
	shufps	$136, %xmm15, %xmm13
	paddd	%xmm1, %xmm13
	movdqa	%xmm13, %xmm1
	psrld	$31, %xmm13
	pand	%xmm0, %xmm1
	paddd	%xmm13, %xmm1
	movdqa	%xmm2, %xmm13
	pand	%xmm0, %xmm13
	psrld	$31, %xmm2
	pandn	%xmm0, %xmm1
	paddd	%xmm13, %xmm2
	paddd	%xmm2, %xmm1
	movdqa	%xmm1, %xmm2
	psrld	$31, %xmm1
	pand	%xmm0, %xmm2
	paddd	%xmm1, %xmm2
	movdqa	%xmm8, %xmm1
	movdqa	%xmm12, %xmm8
	pmuludq	%xmm4, %xmm1
	psrlq	$31, %xmm8
	movdqa	%xmm1, %xmm13
	shufps	$136, %xmm1, %xmm12
	pand	%xmm0, %xmm12
	movdqa	%xmm10, %xmm1
	psrlq	$31, %xmm13
	shufps	$136, %xmm3, %xmm10
	psrlq	$31, %xmm1
	pand	%xmm0, %xmm10
	shufps	$136, %xmm13, %xmm8
	paddd	%xmm12, %xmm8
	movdqa	%xmm3, %xmm12
	psrlq	$31, %xmm12
	movdqa	%xmm8, %xmm3
	pand	%xmm0, %xmm3
	psrld	$31, %xmm8
	shufps	$136, %xmm12, %xmm1
	paddd	%xmm10, %xmm1
	paddd	%xmm8, %xmm3
	movdqa	%xmm1, %xmm8
	pand	%xmm0, %xmm8
	psrld	$31, %xmm1
	paddd	%xmm8, %xmm1
	paddd	%xmm1, %xmm3
	movdqa	%xmm3, %xmm1
	psrld	$31, %xmm3
	pand	%xmm0, %xmm1
	paddd	%xmm3, %xmm1
	movdqa	%xmm2, %xmm3
	pandn	%xmm0, %xmm2
	paddd	%xmm11, %xmm3
	paddd	%xmm11, %xmm2
	movdqa	%xmm3, %xmm8
	psrld	$31, %xmm3
	pand	%xmm0, %xmm8
	paddd	%xmm8, %xmm3
	movups	%xmm3, (%r14,%rsi)
	movdqa	%xmm1, %xmm3
	pandn	%xmm0, %xmm1
	paddd	%xmm9, %xmm3
	paddd	%xmm9, %xmm1
	movdqa	%xmm3, %xmm8
	psrld	$31, %xmm3
	pand	%xmm0, %xmm8
	paddd	%xmm8, %xmm3
	movups	%xmm3, (%rax,%rsi)
	movdqa	%xmm2, %xmm3
	psrld	$31, %xmm2
	pand	%xmm0, %xmm3
	paddd	%xmm3, %xmm2
	movups	%xmm2, (%r8,%rsi)
	movdqa	%xmm1, %xmm2
	psrld	$31, %xmm1
	pand	%xmm0, %xmm2
	paddd	%xmm2, %xmm1
	movups	%xmm1, (%rdi,%rsi)
	addq	$16, %rsi
	cmpq	%rsi, %r9
	jne	.L44
	movl	8(%rsp), %esi
	movq	%rbp, -88(%rsp)
	testl	%esi, %esi
	je	.L46
	movl	16(%rsp), %ebx
	movl	-96(%rsp), %edi
	movq	-64(%rsp), %r9
	movl	-56(%rsp), %r8d
	leal	(%rbx,%rdi), %esi
	movslq	%esi, %rdi
	leaq	(%r9,%rdi,4), %rbx
	movl	(%rbx), %r10d
	movl	%r10d, -16(%rsp)
	movl	-80(%rsp), %r10d
	leal	(%rsi,%r10), %edi
	movslq	%edi, %rdi
	leaq	(%r9,%rdi,4), %r11
	leal	(%rsi,%r8), %edi
	movslq	%edi, %rdi
	leaq	(%r9,%rdi,4), %r10
	movl	(%r10), %r12d
	movl	%r12d, -4(%rsp)
	movl	-76(%rsp), %r12d
	addl	%r12d, %esi
	movl	(%r11), %r12d
	movslq	%esi, %rsi
	leaq	(%r9,%rsi,4), %r9
	movq	%r12, %rsi
	imulq	%rdx, %rsi
	imulq	%rcx, %r12
	movq	%rsi, %rdi
	andl	$2147483647, %esi
	shrq	$31, %rdi
	addl	%edi, %esi
	movl	(%r9), %edi
	movq	%rdi, %r8
	imulq	%rdx, %rdi
	imulq	%rcx, %r8
	movq	%r8, %rbp
	andl	$2147483647, %r8d
	shrq	$31, %rbp
	addl	%ebp, %r8d
	movl	%r8d, %ebp
	shrl	$31, %r8d
	andl	$2147483647, %ebp
	addl	%ebp, %r8d
	movl	%esi, %ebp
	andl	$2147483647, %esi
	shrl	$31, %ebp
	notl	%r8d
	addl	%ebp, %esi
	andl	$2147483647, %r8d
	movl	%r12d, %ebp
	shrq	$31, %r12
	addl	%esi, %r8d
	andl	$2147483647, %ebp
	movl	%r8d, %esi
	shrl	$31, %r8d
	andl	$2147483647, %esi
	addl	%r8d, %esi
	movl	%ebp, %r8d
	movq	%rdi, %rbp
	andl	$2147483647, %edi
	shrq	$31, %rbp
	addl	%r12d, %r8d
	addl	%ebp, %edi
	movl	%r8d, %ebp
	shrl	$31, %r8d
	movl	%edi, %r12d
	andl	$2147483647, %ebp
	andl	$2147483647, %r12d
	addl	%ebp, %r12d
	addl	%r12d, %r8d
	shrl	$31, %edi
	addl	%edi, %r8d
	movl	%r8d, %ebp
	andl	$2147483647, %r8d
	shrl	$31, %ebp
	movl	%ebp, %edi
	addl	%r8d, %edi
	movl	-16(%rsp), %r8d
	addl	%esi, %r8d
	notl	%esi
	movl	%r8d, %ebp
	andl	$2147483647, %r8d
	andl	$2147483647, %esi
	shrl	$31, %ebp
	movl	%ebp, %r12d
	addl	%r8d, %r12d
	movl	%r12d, (%rbx)
	movl	-4(%rsp), %r12d
	leal	(%r12,%rdi), %r8d
	notl	%edi
	movl	%r8d, %ebp
	shrl	$31, %r8d
	andl	$2147483647, %edi
	andl	$2147483647, %ebp
	movl	%ebp, %ebx
	addl	%r8d, %ebx
	movl	%ebx, (%r11)
	movl	-16(%rsp), %r11d
	addl	%r11d, %esi
	movl	%esi, %r11d
	shrl	$31, %esi
	andl	$2147483647, %r11d
	movl	%r11d, %r8d
	addl	%esi, %r8d
	movl	%r12d, %esi
	addl	%edi, %esi
	movl	%r8d, (%r10)
	movl	-80(%rsp), %r10d
	movl	%esi, %r11d
	andl	$2147483647, %esi
	shrl	$31, %r11d
	movl	%r11d, %edi
	addl	%esi, %edi
	movl	12(%rsp), %esi
	movl	%edi, (%r9)
	cmpl	%r10d, %esi
	jge	.L46
	movl	-96(%rsp), %edi
	movq	-64(%rsp), %rbx
	movl	-56(%rsp), %r8d
	addl	%edi, %esi
	movslq	%esi, %rdi
	leaq	(%rbx,%rdi,4), %r11
	leal	(%rsi,%r10), %edi
	movslq	%edi, %rdi
	movl	(%r11), %ebp
	leaq	(%rbx,%rdi,4), %r10
	leal	(%rsi,%r8), %edi
	movslq	%edi, %rdi
	leaq	(%rbx,%rdi,4), %r9
	movl	(%r10), %edi
	movl	(%r9), %r12d
	movl	%r12d, -4(%rsp)
	movl	-76(%rsp), %r12d
	addl	%r12d, %esi
	movslq	%esi, %rsi
	leaq	(%rbx,%rsi,4), %rsi
	movq	%rsi, %rbx
	movq	%rdi, %rsi
	imulq	%rdx, %rsi
	movq	%rbx, -16(%rsp)
	imulq	%rcx, %rdi
	movq	%rsi, %r8
	andl	$2147483647, %esi
	shrq	$31, %r8
	addl	%r8d, %esi
	movl	(%rbx), %r8d
	movq	%rcx, %rbx
	imulq	%r8, %rbx
	imulq	%rdx, %r8
	movq	%rbx, %r12
	andl	$2147483647, %ebx
	shrq	$31, %r12
	addl	%r12d, %ebx
	movl	%ebx, %r12d
	andl	$2147483647, %ebx
	shrl	$31, %r12d
	addl	%r12d, %ebx
	movl	%esi, %r12d
	andl	$2147483647, %esi
	shrl	$31, %r12d
	notl	%ebx
	addl	%r12d, %esi
	andl	$2147483647, %ebx
	addl	%esi, %ebx
	movl	%ebx, %esi
	andl	$2147483647, %ebx
	shrl	$31, %esi
	addl	%ebx, %esi
	movq	%rdi, %rbx
	andl	$2147483647, %edi
	shrq	$31, %rbx
	addl	%ebx, %edi
	movl	%r8d, %ebx
	shrq	$31, %r8
	andl	$2147483647, %ebx
	addl	%ebx, %r8d
	movl	%edi, %ebx
	movl	%r8d, %r12d
	andl	$2147483647, %ebx
	andl	$2147483647, %r12d
	addl	%r12d, %ebx
	shrl	$31, %edi
	addl	%ebx, %edi
	shrl	$31, %r8d
	addl	%edi, %r8d
	movl	%r8d, %edi
	andl	$2147483647, %r8d
	shrl	$31, %edi
	addl	%r8d, %edi
	leal	0(%rbp,%rsi), %r8d
	notl	%esi
	movl	%r8d, %ebx
	andl	$2147483647, %r8d
	andl	$2147483647, %esi
	shrl	$31, %ebx
	addl	%ebp, %esi
	addl	%r8d, %ebx
	movl	%ebx, (%r11)
	movl	-4(%rsp), %ebx
	leal	(%rdi,%rbx), %r8d
	notl	%edi
	movl	%r8d, %r11d
	andl	$2147483647, %r8d
	andl	$2147483647, %edi
	shrl	$31, %r11d
	addl	%r8d, %r11d
	movl	%r11d, (%r10)
	movl	%esi, %r11d
	andl	$2147483647, %esi
	movl	-80(%rsp), %r10d
	shrl	$31, %r11d
	movl	%r11d, %r8d
	addl	%esi, %r8d
	movl	%ebx, %esi
	addl	%edi, %esi
	movl	%r8d, (%r9)
	movl	%esi, %r11d
	andl	$2147483647, %esi
	shrl	$31, %r11d
	movl	%r11d, %edi
	addl	%esi, %edi
	movq	-16(%rsp), %rsi
	movl	%edi, (%rsi)
	movl	20(%rsp), %esi
	cmpl	%r10d, %esi
	jge	.L46
	movl	-96(%rsp), %edi
	movq	-64(%rsp), %r9
	movl	-56(%rsp), %r8d
	movl	-76(%rsp), %r12d
	addl	%edi, %esi
	movslq	%esi, %rdi
	leaq	(%r9,%rdi,4), %rbx
	leal	(%r10,%rsi), %edi
	movslq	%edi, %rdi
	movl	(%rbx), %ebp
	leaq	(%r9,%rdi,4), %r11
	leal	(%r8,%rsi), %edi
	addl	%r12d, %esi
	movslq	%edi, %rdi
	movslq	%esi, %rsi
	leaq	(%r9,%rdi,4), %r10
	leaq	(%r9,%rsi,4), %r8
	movl	(%r11), %esi
	movl	(%r10), %edi
	movq	%r8, -16(%rsp)
	movl	(%r8), %r8d
	movl	%edi, -4(%rsp)
	movq	%rsi, %rdi
	imulq	%rcx, %rsi
	imulq	%rdx, %rdi
	movl	%edi, %r9d
	shrq	$31, %rdi
	andl	$2147483647, %r9d
	addl	%edi, %r9d
	movq	%r8, %rdi
	imulq	%rcx, %rdi
	imulq	%rdx, %r8
	movl	%edi, %r12d
	shrq	$31, %rdi
	andl	$2147483647, %r12d
	movl	%r8d, %edx
	shrq	$31, %r8
	addl	%r12d, %edi
	andl	$2147483647, %edx
	movl	%edi, %r12d
	shrl	$31, %edi
	addl	%r8d, %edx
	andl	$2147483647, %r12d
	addl	%r12d, %edi
	movl	%r9d, %r12d
	shrl	$31, %r9d
	andl	$2147483647, %r12d
	notl	%edi
	addl	%r12d, %r9d
	andl	$2147483647, %edi
	addl	%r9d, %edi
	movl	%edi, %r9d
	shrl	$31, %edi
	andl	$2147483647, %r9d
	addl	%r9d, %edi
	movl	%esi, %r9d
	shrq	$31, %rsi
	andl	$2147483647, %r9d
	leal	(%r9,%rsi), %ecx
	movl	%edx, %esi
	andl	$2147483647, %esi
	movl	%esi, %r8d
	movl	%ecx, %esi
	shrl	$31, %ecx
	andl	$2147483647, %esi
	addl	%r8d, %esi
	movq	-16(%rsp), %r8
	addl	%esi, %ecx
	shrl	$31, %edx
	addl	%edx, %ecx
	movl	%ecx, %edx
	shrl	$31, %ecx
	andl	$2147483647, %edx
	addl	%ecx, %edx
	leal	(%rdi,%rbp), %ecx
	notl	%edi
	movl	%ecx, %esi
	shrl	$31, %ecx
	andl	$2147483647, %edi
	andl	$2147483647, %esi
	addl	%ecx, %esi
	movl	%esi, (%rbx)
	movl	-4(%rsp), %ebx
	leal	(%rdx,%rbx), %ecx
	notl	%edx
	movl	%ecx, %esi
	shrl	$31, %ecx
	andl	$2147483647, %edx
	andl	$2147483647, %esi
	addl	%ebx, %edx
	addl	%ecx, %esi
	leal	(%rdi,%rbp), %ecx
	movl	%esi, (%r11)
	movl	%ecx, %esi
	shrl	$31, %ecx
	andl	$2147483647, %esi
	addl	%ecx, %esi
	movl	%edx, %ecx
	shrl	$31, %edx
	andl	$2147483647, %ecx
	movl	%esi, (%r10)
	addl	%edx, %ecx
	movl	%ecx, (%r8)
.L46:
	movq	-24(%rsp), %rcx
	addq	$8, -104(%rsp)
	addq	%rcx, -72(%rsp)
	movq	-104(%rsp), %rdi
	addq	%rcx, -112(%rsp)
	addq	%rcx, %r14
	addq	%rcx, %rax
	addq	%rcx, -120(%rsp)
	cmpq	%rdi, -32(%rsp)
	jne	.L47
	movq	-64(%rsp), %r12
.L48:
	sarl	-80(%rsp)
	sall	32(%rsp)
	subl	$1, 36(%rsp)
	jne	.L38
	addq	$48, %rsp
	.cfi_remember_state
	.cfi_def_cfa_offset 56
	popq	%rbx
	.cfi_def_cfa_offset 48
	popq	%rbp
	.cfi_def_cfa_offset 40
	popq	%r12
	.cfi_def_cfa_offset 32
	popq	%r13
	.cfi_def_cfa_offset 24
	popq	%r14
	.cfi_def_cfa_offset 16
	popq	%r15
	.cfi_def_cfa_offset 8
	ret
	.p2align 4,,10
	.p2align 3
.L54:
	.cfi_restore_state
	movq	-88(%rsp), %rbp
	movq	%r14, %rdi
	.p2align 4,,10
	.p2align 3
.L43:
	movl	(%rdi,%rbp,4), %r11d
	movl	(%rdi), %r10d
	movl	(%rdi,%r15,4), %r9d
	movq	%r11, %rsi
	imulq	%rcx, %r11
	imulq	%rdx, %rsi
	movl	%esi, %r8d
	shrq	$31, %rsi
	andl	$2147483647, %r8d
	addl	%esi, %r8d
	movl	(%rdi,%r13,4), %esi
	movq	%rsi, %r12
	imulq	%rdx, %rsi
	imulq	%rcx, %r12
	movl	%r12d, %ebx
	shrq	$31, %r12
	andl	$2147483647, %ebx
	addl	%r12d, %ebx
	movl	%ebx, %r12d
	shrl	$31, %ebx
	andl	$2147483647, %r12d
	addl	%r12d, %ebx
	movl	%r8d, %r12d
	shrl	$31, %r8d
	andl	$2147483647, %r12d
	notl	%ebx
	addl	%r12d, %r8d
	andl	$2147483647, %ebx
	addl	%r8d, %ebx
	movl	%ebx, %r8d
	shrl	$31, %ebx
	andl	$2147483647, %r8d
	addl	%ebx, %r8d
	movl	%r11d, %ebx
	shrq	$31, %r11
	andl	$2147483647, %ebx
	addl	%ebx, %r11d
	movl	%esi, %ebx
	shrq	$31, %rsi
	andl	$2147483647, %ebx
	movl	%r11d, %r12d
	shrl	$31, %r11d
	addl	%ebx, %esi
	andl	$2147483647, %r12d
	movl	%esi, %ebx
	andl	$2147483647, %ebx
	addl	%r12d, %ebx
	addl	%ebx, %r11d
	shrl	$31, %esi
	addl	%esi, %r11d
	movl	%r11d, %esi
	shrl	$31, %r11d
	andl	$2147483647, %esi
	addl	%r11d, %esi
	leal	(%r8,%r10), %r11d
	notl	%r8d
	movl	%r11d, %ebx
	shrl	$31, %r11d
	andl	$2147483647, %r8d
	andl	$2147483647, %ebx
	addl	%r10d, %r8d
	addl	%ebx, %r11d
	movl	%r8d, %r10d
	shrl	$31, %r8d
	movl	%r11d, (%rdi)
	leal	(%rsi,%r9), %r11d
	notl	%esi
	andl	$2147483647, %r10d
	movl	%r11d, %ebx
	andl	$2147483647, %esi
	shrl	$31, %r11d
	addl	%r10d, %r8d
	andl	$2147483647, %ebx
	addl	%r9d, %esi
	addl	%ebx, %r11d
	movl	%r11d, (%rdi,%rbp,4)
	movl	%r8d, (%rdi,%r15,4)
	movl	%esi, %r8d
	shrl	$31, %esi
	andl	$2147483647, %r8d
	addl	%r8d, %esi
	movl	%esi, (%rdi,%r13,4)
	addq	$4, %rdi
	cmpq	%rax, %rdi
	jne	.L43
	movq	%rbp, -88(%rsp)
	jmp	.L46
	.p2align 4,,10
	.p2align 3
.L3:
	movslq	-76(%rsp), %rdi
	movq	%r12, %rcx
	leaq	(%r12,%rbp,4), %r11
	.p2align 4,,10
	.p2align 3
.L9:
	movl	(%rcx), %r9d
	movl	(%rcx,%rbp,4), %edx
	movl	(%rcx,%r15,4), %esi
	movl	(%rcx,%rdi,4), %eax
	leal	(%r9,%rdx), %r8d
	notl	%edx
	movl	%r8d, %r10d
	andl	$2147483647, %r8d
	andl	$2147483647, %edx
	shrl	$31, %r10d
	addl	%r9d, %edx
	addl	%r10d, %r8d
	movl	%r8d, (%rcx)
	leal	(%rsi,%rax), %r8d
	notl	%eax
	movl	%r8d, %r10d
	andl	$2147483647, %r8d
	andl	$2147483647, %eax
	shrl	$31, %r10d
	addl	%esi, %eax
	addl	%r10d, %r8d
	movl	%r8d, (%rcx,%rbp,4)
	movl	%edx, %r8d
	andl	$2147483647, %edx
	shrl	$31, %r8d
	addl	%r8d, %edx
	movl	%edx, (%rcx,%r15,4)
	movl	%eax, %edx
	andl	$2147483647, %eax
	shrl	$31, %edx
	addl	%edx, %eax
	movl	%eax, (%rcx,%rdi,4)
	addq	$4, %rcx
	cmpq	%r11, %rcx
	jne	.L9
	cmpl	$1, -80(%rsp)
	movslq	-52(%rsp), %r14
	jne	.L6
.L11:
	leaq	(%r14,%rbp), %rax
	leaq	(%r12,%r14,4), %rdx
	leaq	(%r12,%rax,4), %r11
	.p2align 4,,10
	.p2align 3
.L18:
	movl	(%rdx,%rdi,4), %r9d
	movl	(%rdx), %r8d
	movl	(%rdx,%rbp,4), %eax
	movl	(%rdx,%r15,4), %esi
	movl	%r9d, %ecx
	notl	%ecx
	andl	$2147483647, %ecx
	addl	%r8d, %ecx
	movl	%ecx, %r10d
	andl	$2147483647, %ecx
	shrl	$31, %r10d
	addl	%r10d, %ecx
	movl	%ecx, (%rdx)
	leal	(%rax,%rsi), %ecx
	notl	%eax
	movl	%ecx, %r10d
	andl	$2147483647, %ecx
	andl	$2147483647, %eax
	shrl	$31, %r10d
	addl	%esi, %eax
	addl	%r10d, %ecx
	movl	%ecx, (%rdx,%rbp,4)
	leal	(%r8,%r9), %ecx
	movl	%ecx, %r8d
	andl	$2147483647, %ecx
	shrl	$31, %r8d
	addl	%r8d, %ecx
	movl	%ecx, (%rdx,%r15,4)
	movl	%eax, %ecx
	andl	$2147483647, %eax
	shrl	$31, %ecx
	addl	%ecx, %eax
	movl	%eax, (%rdx,%rdi,4)
	addq	$4, %rdx
	cmpq	%r11, %rdx
	jne	.L18
	cmpl	$1, -80(%rsp)
	movl	$8, -120(%rsp)
	jne	.L138
.L20:
	movslq	-120(%rsp), %rax
	leaq	(%r12,%rax,4), %rcx
	addq	%rbp, %rax
	leaq	(%r12,%rax,4), %r9
	.p2align 4,,10
	.p2align 3
.L27:
	movl	(%rcx,%rbp,4), %eax
	movl	(%rcx,%rdi,4), %edx
	movl	(%rcx), %r10d
	movl	(%rcx,%r15,4), %esi
	movl	%eax, %r8d
	shrl	$16, %eax
	sall	$15, %r8d
	andl	$2147450880, %r8d
	orl	%eax, %r8d
	movl	%edx, %eax
	shrl	$16, %edx
	sall	$15, %eax
	andl	$2147450880, %eax
	orl	%edx, %eax
	movl	%eax, %r11d
	xorl	$2147483647, %r11d
	addl	%r8d, %r11d
	addl	%eax, %r8d
	movl	%r11d, %edx
	movl	%r8d, %eax
	andl	$2147483647, %r11d
	andl	$2147483647, %r8d
	shrl	$31, %edx
	shrl	$31, %eax
	addl	%r11d, %edx
	addl	%r8d, %eax
	leal	(%r10,%rdx), %r8d
	notl	%edx
	movl	%r8d, %r11d
	andl	$2147483647, %r8d
	andl	$2147483647, %edx
	shrl	$31, %r11d
	addl	%r10d, %edx
	addl	%r11d, %r8d
	movl	%r8d, (%rcx)
	leal	(%rsi,%rax), %r8d
	notl	%eax
	movl	%r8d, %r11d
	andl	$2147483647, %r8d
	andl	$2147483647, %eax
	shrl	$31, %r11d
	addl	%esi, %eax
	addl	%r11d, %r8d
	movl	%r8d, (%rcx,%rbp,4)
	movl	%edx, %r8d
	andl	$2147483647, %edx
	shrl	$31, %r8d
	addl	%r8d, %edx
	movl	%edx, (%rcx,%r15,4)
	movl	%eax, %edx
	andl	$2147483647, %eax
	shrl	$31, %edx
	addl	%edx, %eax
	movl	%eax, (%rcx,%rdi,4)
	addq	$4, %rcx
	cmpq	%rcx, %r9
	jne	.L27
	movl	-120(%rsp), %eax
	movl	-52(%rsp), %ecx
	addl	%ecx, %eax
	cmpl	$1, -80(%rsp)
	movl	%eax, -96(%rsp)
	jne	.L24
.L29:
	movslq	-96(%rsp), %rax
	leaq	(%r12,%rax,4), %rdx
	addq	%rbp, %rax
	leaq	(%r12,%rax,4), %r10
	.p2align 4,,10
	.p2align 3
.L35:
	movl	(%rdx,%rbp,4), %eax
	movl	(%rdx,%rdi,4), %r8d
	movl	(%rdx), %esi
	movl	(%rdx,%r15,4), %ecx
	movl	%eax, %r9d
	shrl	$16, %eax
	sall	$15, %r9d
	andl	$2147450880, %r9d
	orl	%eax, %r9d
	movl	%r8d, %eax
	shrl	$16, %r8d
	sall	$15, %eax
	andl	$2147450880, %eax
	orl	%r8d, %eax
	leal	(%rax,%r9), %r11d
	xorl	$2147483647, %eax
	addl	%r9d, %eax
	movl	%r11d, %r8d
	shrl	$31, %r11d
	andl	$2147483647, %r8d
	movl	%eax, %r9d
	shrl	$31, %eax
	addl	%r11d, %r8d
	andl	$2147483647, %r9d
	addl	%r9d, %eax
	movl	%r8d, %r9d
	addl	%esi, %r8d
	notl	%r9d
	andl	$2147483647, %r9d
	addl	%esi, %r9d
	movl	%r8d, %esi
	shrl	$31, %r8d
	movl	%r9d, %r11d
	shrl	$31, %r9d
	andl	$2147483647, %esi
	andl	$2147483647, %r11d
	addl	%r8d, %esi
	addl	%r11d, %r9d
	movl	%r9d, (%rdx)
	leal	(%rax,%rcx), %r9d
	notl	%eax
	andl	$2147483647, %eax
	movl	%r9d, %r11d
	shrl	$31, %r9d
	addl	%ecx, %eax
	andl	$2147483647, %r11d
	movl	%eax, %ecx
	shrl	$31, %eax
	addl	%r11d, %r9d
	andl	$2147483647, %ecx
	movl	%r9d, (%rdx,%rbp,4)
	addl	%ecx, %eax
	movl	%esi, (%rdx,%r15,4)
	movl	%eax, (%rdx,%rdi,4)
	addq	$4, %rdx
	cmpq	%rdx, %r10
	jne	.L35
	jmp	.L36
.L139:
	movslq	-76(%rsp), %rdi
	jmp	.L20
.L140:
	movslq	-76(%rsp), %rdi
	jmp	.L29
.L137:
	movslq	-76(%rsp), %rdi
	jmp	.L11
.L49:
	xorl	%edx, %edx
	xorl	%eax, %eax
	jmp	.L4
.L52:
	xorl	%edx, %edx
	xorl	%eax, %eax
	jmp	.L21
.L53:
	xorl	%edx, %edx
	xorl	%eax, %eax
	jmp	.L30
.L50:
	xorl	%eax, %eax
	xorl	%edx, %edx
	jmp	.L12
	.cfi_endproc
.LFE49:
	.size	mrsn_ntt_256, .-mrsn_ntt_256
	.p2align 4
	.globl	mrsn_invntt_256
	.type	mrsn_invntt_256, @function
mrsn_invntt_256:
.LFB50:
	.cfi_startproc
	endbr64
	pushq	%r15
	.cfi_def_cfa_offset 16
	.cfi_offset 15, -16
	pushq	%r14
	.cfi_def_cfa_offset 24
	.cfi_offset 14, -24
	movq	%rdi, %r14
	pushq	%r13
	.cfi_def_cfa_offset 32
	.cfi_offset 13, -32
	pushq	%r12
	.cfi_def_cfa_offset 40
	.cfi_offset 12, -40
	pushq	%rbp
	.cfi_def_cfa_offset 48
	.cfi_offset 6, -48
	pushq	%rbx
	.cfi_def_cfa_offset 56
	.cfi_offset 3, -56
	subq	$64, %rsp
	.cfi_def_cfa_offset 120
	movdqa	.LC0(%rip), %xmm0
	movq	.LC1(%rip), %xmm14
	movq	%rdi, 48(%rsp)
	movdqa	.LC3(%rip), %xmm15
	movl	$5, 44(%rsp)
	movl	$64, 40(%rsp)
	movl	$1, -40(%rsp)
	.p2align 4,,10
	.p2align 3
.L183:
	movl	-40(%rsp), %esi
	movl	%esi, %eax
	movl	%esi, -72(%rsp)
	addl	%esi, %esi
	movl	%esi, -40(%rsp)
	addl	%eax, %esi
	movl	%esi, -68(%rsp)
	leal	0(,%rax,4), %esi
	movl	%esi, -36(%rsp)
	leal	0(,%rax,8), %esi
	movl	%esi, -120(%rsp)
	testl	%eax, %eax
	jle	.L143
	movslq	-40(%rsp), %r15
	movl	$1, %r13d
	cmpl	$1, %eax
	je	.L144
	movslq	-72(%rsp), %r13
	leaq	0(,%r15,4), %rcx
	leaq	16(%r14), %r10
	leaq	16(%rcx), %rbp
	leaq	(%r14,%rcx), %r9
	leaq	0(,%r13,4), %rdx
	leaq	(%r15,%r13), %r12
	leaq	16(%rdx), %r11
	cmpq	%rdx, %rbp
	leaq	(%r14,%rdx), %r8
	setle	%dil
	leaq	16(,%r12,4), %rax
	cmpq	%r11, %rcx
	movq	%rax, -104(%rsp)
	setge	%al
	leaq	0(,%r12,4), %rsi
	orl	%edi, %eax
	cmpq	%r10, %r9
	leaq	(%r14,%rsi), %rbx
	setnb	%dil
	andl	%eax, %edi
	cmpq	%r10, %r8
	setnb	%al
	andl	%edi, %eax
	cmpq	%rcx, -104(%rsp)
	setle	%cl
	cmpq	%rsi, %rbp
	setle	%dil
	orl	%edi, %ecx
	andl	%ecx, %eax
	cmpq	%rdx, -104(%rsp)
	setle	%dl
	cmpq	%r11, %rsi
	setge	%cl
	orl	%ecx, %edx
	testb	%dl, %al
	je	.L144
	cmpq	%r10, %rbx
	jb	.L144
	movl	-72(%rsp), %esi
	leal	-1(%rsi), %eax
	movl	%esi, %ecx
	cmpl	$2, %eax
	jbe	.L194
	movdqu	(%r14), %xmm6
	movdqu	(%r9), %xmm2
	movl	%esi, %eax
	movdqu	(%r8), %xmm4
	movdqu	(%rbx), %xmm1
	shrl	$2, %eax
	movdqa	%xmm6, %xmm3
	paddd	%xmm2, %xmm3
	pandn	%xmm0, %xmm2
	movdqa	%xmm3, %xmm5
	pand	%xmm0, %xmm3
	paddd	%xmm6, %xmm2
	psrld	$31, %xmm5
	paddd	%xmm5, %xmm3
	movups	%xmm3, (%r14)
	movdqa	%xmm4, %xmm3
	paddd	%xmm1, %xmm3
	pandn	%xmm0, %xmm1
	movdqa	%xmm3, %xmm5
	pand	%xmm0, %xmm3
	paddd	%xmm4, %xmm1
	psrld	$31, %xmm5
	paddd	%xmm5, %xmm3
	movups	%xmm3, (%r9)
	movdqa	%xmm2, %xmm3
	pand	%xmm0, %xmm2
	psrld	$31, %xmm3
	paddd	%xmm3, %xmm2
	movups	%xmm2, (%r8)
	movdqa	%xmm1, %xmm2
	pand	%xmm0, %xmm1
	psrld	$31, %xmm2
	paddd	%xmm2, %xmm1
	movups	%xmm1, (%rbx)
	cmpl	$1, %eax
	je	.L146
	movdqu	16(%r14), %xmm6
	movdqu	16(%r9), %xmm2
	movdqu	16(%r8), %xmm4
	movdqu	16(%rbx), %xmm1
	movdqa	%xmm6, %xmm3
	paddd	%xmm2, %xmm3
	pandn	%xmm0, %xmm2
	movdqa	%xmm3, %xmm5
	pand	%xmm0, %xmm3
	paddd	%xmm6, %xmm2
	psrld	$31, %xmm5
	paddd	%xmm5, %xmm3
	movups	%xmm3, 16(%r14)
	movdqa	%xmm4, %xmm3
	paddd	%xmm1, %xmm3
	pandn	%xmm0, %xmm1
	movdqa	%xmm3, %xmm5
	pand	%xmm0, %xmm3
	paddd	%xmm4, %xmm1
	psrld	$31, %xmm5
	paddd	%xmm5, %xmm3
	movups	%xmm3, 16(%r9)
	movdqa	%xmm2, %xmm3
	pand	%xmm0, %xmm2
	psrld	$31, %xmm3
	paddd	%xmm3, %xmm2
	movups	%xmm2, 16(%r8)
	movdqa	%xmm1, %xmm2
	pand	%xmm0, %xmm1
	psrld	$31, %xmm2
	paddd	%xmm2, %xmm1
	movups	%xmm1, 16(%rbx)
	cmpl	$2, %eax
	je	.L146
	movdqu	32(%r14), %xmm6
	movdqu	32(%r9), %xmm2
	movdqu	32(%r8), %xmm4
	movdqu	32(%rbx), %xmm1
	movdqa	%xmm6, %xmm3
	paddd	%xmm2, %xmm3
	pandn	%xmm0, %xmm2
	movdqa	%xmm3, %xmm5
	pand	%xmm0, %xmm3
	paddd	%xmm6, %xmm2
	psrld	$31, %xmm5
	paddd	%xmm5, %xmm3
	movups	%xmm3, 32(%r14)
	movdqa	%xmm4, %xmm3
	paddd	%xmm1, %xmm3
	pandn	%xmm0, %xmm1
	movdqa	%xmm3, %xmm5
	pand	%xmm0, %xmm3
	paddd	%xmm4, %xmm1
	psrld	$31, %xmm5
	paddd	%xmm5, %xmm3
	movups	%xmm3, 32(%r9)
	movdqa	%xmm2, %xmm3
	pand	%xmm0, %xmm2
	psrld	$31, %xmm3
	paddd	%xmm3, %xmm2
	movups	%xmm2, 32(%r8)
	movdqa	%xmm1, %xmm2
	pand	%xmm0, %xmm1
	psrld	$31, %xmm2
	paddd	%xmm2, %xmm1
	movups	%xmm1, 32(%rbx)
	cmpl	$4, %eax
	jne	.L146
	movdqu	48(%r14), %xmm6
	movdqu	48(%r9), %xmm2
	movdqu	48(%r8), %xmm4
	movdqu	48(%rbx), %xmm1
	movdqa	%xmm6, %xmm3
	paddd	%xmm2, %xmm3
	pandn	%xmm0, %xmm2
	movdqa	%xmm3, %xmm5
	pand	%xmm0, %xmm3
	paddd	%xmm6, %xmm2
	psrld	$31, %xmm5
	paddd	%xmm5, %xmm3
	movups	%xmm3, 48(%r14)
	movdqa	%xmm4, %xmm3
	paddd	%xmm1, %xmm3
	pandn	%xmm0, %xmm1
	movdqa	%xmm3, %xmm5
	pand	%xmm0, %xmm3
	paddd	%xmm4, %xmm1
	psrld	$31, %xmm5
	paddd	%xmm5, %xmm3
	movups	%xmm3, 48(%r9)
	movdqa	%xmm2, %xmm3
	pand	%xmm0, %xmm2
	psrld	$31, %xmm3
	paddd	%xmm3, %xmm2
	movups	%xmm2, 48(%r8)
	movdqa	%xmm1, %xmm2
	pand	%xmm0, %xmm1
	psrld	$31, %xmm2
	paddd	%xmm2, %xmm1
	movups	%xmm1, 48(%rbx)
	.p2align 4,,10
	.p2align 3
.L149:
	movslq	-36(%rsp), %r12
	leaq	(%r15,%r12), %rbp
	leaq	0(%r13,%r12), %rdx
	leaq	0(,%rbp,4), %rsi
	leaq	0(,%rdx,4), %rdi
	movq	%rdx, -64(%rsp)
	leaq	16(%rsi), %r10
	leaq	16(%rdi), %r9
	cmpq	%rdi, %r10
	leaq	0(,%r12,4), %r8
	leaq	0(%rbp,%r13), %rbx
	setle	%dl
	cmpq	%r9, %rsi
	leaq	16(%r8), %r11
	setge	%cl
	leaq	16(,%rbx,4), %rax
	orl	%ecx, %edx
	cmpq	%r8, %r10
	movq	%rax, -104(%rsp)
	leaq	0(,%rbx,4), %rax
	setle	%cl
	cmpq	%r11, %rsi
	setge	-80(%rsp)
	orb	-80(%rsp), %cl
	andl	%edx, %ecx
	cmpq	%r8, %r9
	setle	%dl
	cmpq	%r11, %rdi
	setge	-80(%rsp)
	orb	-80(%rsp), %dl
	andl	%ecx, %edx
	cmpq	%rsi, -104(%rsp)
	setle	%cl
	cmpq	%rax, %r10
	setle	%r10b
	orl	%r10d, %ecx
	andl	%edx, %ecx
	cmpq	%rdi, -104(%rsp)
	setle	%dl
	cmpq	%r9, %rax
	setge	%r9b
	orl	%r9d, %edx
	testb	%dl, %cl
	je	.L302
	cmpq	%r8, -104(%rsp)
	setle	%dl
	cmpq	%r11, %rax
	setge	%cl
	orb	%dl, %cl
	je	.L302
	movl	-72(%rsp), %r9d
	leal	-1(%r9), %edx
	movl	%r9d, %ecx
	cmpl	$2, %edx
	jbe	.L195
	leaq	(%r14,%r8), %rdx
	leaq	(%r14,%rsi), %rcx
	addq	%r14, %rdi
	addq	%r14, %rax
	movdqu	(%rdx), %xmm1
	movdqu	(%rcx), %xmm4
	shrl	$2, %r9d
	movdqu	(%rdi), %xmm6
	movdqu	(%rax), %xmm2
	movdqa	%xmm1, %xmm3
	pandn	%xmm0, %xmm1
	paddd	%xmm4, %xmm3
	paddd	%xmm4, %xmm1
	movdqa	%xmm3, %xmm5
	pand	%xmm0, %xmm3
	psrld	$31, %xmm5
	paddd	%xmm5, %xmm3
	movups	%xmm3, (%rdx)
	movdqa	%xmm6, %xmm3
	paddd	%xmm2, %xmm3
	pandn	%xmm0, %xmm2
	movdqa	%xmm3, %xmm5
	pand	%xmm0, %xmm3
	paddd	%xmm6, %xmm2
	psrld	$31, %xmm5
	paddd	%xmm5, %xmm3
	movups	%xmm3, (%rcx)
	movdqa	%xmm2, %xmm3
	pand	%xmm0, %xmm2
	psrld	$31, %xmm3
	paddd	%xmm3, %xmm2
	movups	%xmm2, (%rdi)
	movdqa	%xmm1, %xmm2
	psrld	$31, %xmm1
	pand	%xmm0, %xmm2
	paddd	%xmm2, %xmm1
	movups	%xmm1, (%rax)
	cmpl	$1, %r9d
	je	.L156
	movdqu	16(%rdx), %xmm1
	movdqu	16(%rcx), %xmm4
	movdqu	16(%rdi), %xmm6
	movdqu	16(%rax), %xmm2
	movdqa	%xmm1, %xmm3
	pandn	%xmm0, %xmm1
	paddd	%xmm4, %xmm3
	paddd	%xmm4, %xmm1
	movdqa	%xmm3, %xmm5
	psrld	$31, %xmm3
	pand	%xmm0, %xmm5
	paddd	%xmm5, %xmm3
	movups	%xmm3, 16(%rdx)
	movdqa	%xmm6, %xmm3
	paddd	%xmm2, %xmm3
	pandn	%xmm0, %xmm2
	movdqa	%xmm3, %xmm5
	pand	%xmm0, %xmm3
	paddd	%xmm6, %xmm2
	psrld	$31, %xmm5
	paddd	%xmm5, %xmm3
	movups	%xmm3, 16(%rcx)
	movdqa	%xmm2, %xmm3
	pand	%xmm0, %xmm2
	psrld	$31, %xmm3
	paddd	%xmm3, %xmm2
	movups	%xmm2, 16(%rdi)
	movdqa	%xmm1, %xmm2
	pand	%xmm0, %xmm1
	psrld	$31, %xmm2
	paddd	%xmm2, %xmm1
	movups	%xmm1, 16(%rax)
	cmpl	$2, %r9d
	je	.L156
	movdqu	32(%rdx), %xmm1
	movdqu	32(%rcx), %xmm4
	movdqu	32(%rdi), %xmm6
	movdqu	32(%rax), %xmm2
	movdqa	%xmm1, %xmm3
	pandn	%xmm0, %xmm1
	paddd	%xmm4, %xmm3
	paddd	%xmm4, %xmm1
	movdqa	%xmm3, %xmm5
	pand	%xmm0, %xmm3
	psrld	$31, %xmm5
	paddd	%xmm5, %xmm3
	movups	%xmm3, 32(%rdx)
	movdqa	%xmm6, %xmm3
	paddd	%xmm2, %xmm3
	pandn	%xmm0, %xmm2
	movdqa	%xmm3, %xmm5
	pand	%xmm0, %xmm3
	paddd	%xmm6, %xmm2
	psrld	$31, %xmm5
	paddd	%xmm5, %xmm3
	movups	%xmm3, 32(%rcx)
	movdqa	%xmm2, %xmm3
	pand	%xmm0, %xmm2
	psrld	$31, %xmm3
	paddd	%xmm3, %xmm2
	movups	%xmm2, 32(%rdi)
	movdqa	%xmm1, %xmm2
	pand	%xmm0, %xmm1
	psrld	$31, %xmm2
	paddd	%xmm2, %xmm1
	movups	%xmm1, 32(%rax)
	cmpl	$4, %r9d
	jne	.L156
	movdqu	48(%rdx), %xmm1
	movdqu	48(%rcx), %xmm4
	movdqu	48(%rdi), %xmm6
	movdqu	48(%rax), %xmm2
	movdqa	%xmm1, %xmm3
	pandn	%xmm0, %xmm1
	paddd	%xmm4, %xmm3
	paddd	%xmm4, %xmm1
	movdqa	%xmm3, %xmm5
	pand	%xmm0, %xmm3
	psrld	$31, %xmm5
	paddd	%xmm5, %xmm3
	movups	%xmm3, 48(%rdx)
	movdqa	%xmm6, %xmm3
	paddd	%xmm2, %xmm3
	pandn	%xmm0, %xmm2
	movdqa	%xmm3, %xmm5
	pand	%xmm0, %xmm3
	paddd	%xmm6, %xmm2
	psrld	$31, %xmm5
	paddd	%xmm5, %xmm3
	movups	%xmm3, 48(%rcx)
	movdqa	%xmm2, %xmm3
	pand	%xmm0, %xmm2
	psrld	$31, %xmm3
	paddd	%xmm3, %xmm2
	movups	%xmm2, 48(%rdi)
	movdqa	%xmm1, %xmm2
	pand	%xmm0, %xmm1
	psrld	$31, %xmm2
	paddd	%xmm2, %xmm1
	movups	%xmm1, 48(%rax)
	.p2align 4,,10
	.p2align 3
.L159:
	leaq	(%r12,%r12), %rbp
	leaq	(%r15,%r13), %rax
	addq	%rbp, %rax
	leaq	(%r15,%rbp), %rsi
	leaq	0(,%r12,8), %rdi
	leaq	16(,%rax,4), %r9
	movq	%rsi, -64(%rsp)
	leaq	16(%rdi), %r11
	leaq	0(,%rsi,4), %rcx
	leaq	0(%r13,%rbp), %rsi
	movq	%rax, -80(%rsp)
	leaq	16(%rcx), %r8
	leaq	0(,%rax,4), %rax
	movq	%rsi, -56(%rsp)
	salq	$2, %rsi
	cmpq	%rcx, %r9
	leaq	16(%rsi), %rbx
	setle	%r10b
	cmpq	%rax, %r8
	setle	%dl
	orl	%edx, %r10d
	cmpq	%rsi, %r9
	setle	%dl
	cmpq	%rbx, %rax
	setge	-104(%rsp)
	orb	-104(%rsp), %dl
	andl	%r10d, %edx
	cmpq	%rdi, %r9
	setle	%r9b
	cmpq	%r11, %rax
	setge	%r10b
	orl	%r10d, %r9d
	andl	%edx, %r9d
	cmpq	%rsi, %r8
	setle	%dl
	cmpq	%rbx, %rcx
	setge	%r10b
	orl	%r10d, %edx
	andl	%r9d, %edx
	cmpq	%rdi, %r8
	setle	%r8b
	cmpq	%r11, %rcx
	setge	%r9b
	orl	%r9d, %r8d
	testb	%r8b, %dl
	je	.L303
	cmpq	%rdi, %rbx
	setle	%dl
	cmpq	%r11, %rsi
	setge	%r8b
	orb	%dl, %r8b
	je	.L303
	movl	-72(%rsp), %ebx
	leal	-1(%rbx), %edx
	movl	%ebx, %r8d
	cmpl	$2, %edx
	jbe	.L196
	leaq	(%r14,%rdi), %rdx
	addq	%r14, %rcx
	addq	%r14, %rsi
	addq	%r14, %rax
	movdqu	(%rdx), %xmm6
	movdqu	(%rcx), %xmm2
	shrl	$2, %ebx
	movdqu	(%rsi), %xmm4
	movdqu	(%rax), %xmm1
	movdqa	%xmm6, %xmm3
	paddd	%xmm2, %xmm3
	pandn	%xmm0, %xmm2
	movdqa	%xmm3, %xmm5
	pand	%xmm0, %xmm3
	paddd	%xmm6, %xmm2
	psrld	$31, %xmm5
	paddd	%xmm5, %xmm3
	movups	%xmm3, (%rdx)
	movdqa	%xmm4, %xmm3
	paddd	%xmm1, %xmm3
	pandn	%xmm0, %xmm1
	movdqa	%xmm3, %xmm5
	pand	%xmm0, %xmm3
	paddd	%xmm4, %xmm1
	psrld	$31, %xmm5
	paddd	%xmm5, %xmm3
	movups	%xmm3, (%rcx)
	movdqa	%xmm2, %xmm3
	psrld	$31, %xmm2
	pand	%xmm0, %xmm3
	paddd	%xmm3, %xmm2
	movdqa	%xmm1, %xmm3
	pand	%xmm0, %xmm3
	psrld	$31, %xmm1
	paddd	%xmm1, %xmm3
	movdqa	%xmm2, %xmm1
	pslld	$15, %xmm1
	psrld	$16, %xmm2
	pand	%xmm15, %xmm1
	por	%xmm2, %xmm1
	movdqa	%xmm3, %xmm2
	pslld	$15, %xmm2
	psrld	$16, %xmm3
	pand	%xmm15, %xmm2
	por	%xmm3, %xmm2
	movdqa	%xmm2, %xmm3
	paddd	%xmm1, %xmm3
	pxor	%xmm0, %xmm1
	paddd	%xmm2, %xmm1
	movdqa	%xmm3, %xmm4
	movdqa	%xmm1, %xmm2
	pand	%xmm0, %xmm4
	pand	%xmm0, %xmm1
	psrld	$31, %xmm3
	psrld	$31, %xmm2
	paddd	%xmm4, %xmm3
	paddd	%xmm2, %xmm1
	movups	%xmm3, (%rsi)
	movups	%xmm1, (%rax)
	cmpl	$1, %ebx
	je	.L166
	movdqu	16(%rdx), %xmm6
	movdqu	16(%rcx), %xmm2
	movdqu	16(%rsi), %xmm4
	movdqu	16(%rax), %xmm1
	movdqa	%xmm6, %xmm3
	paddd	%xmm2, %xmm3
	pandn	%xmm0, %xmm2
	movdqa	%xmm3, %xmm5
	pand	%xmm0, %xmm3
	paddd	%xmm6, %xmm2
	psrld	$31, %xmm5
	paddd	%xmm5, %xmm3
	movups	%xmm3, 16(%rdx)
	movdqa	%xmm4, %xmm3
	paddd	%xmm1, %xmm3
	pandn	%xmm0, %xmm1
	movdqa	%xmm3, %xmm5
	psrld	$31, %xmm3
	paddd	%xmm4, %xmm1
	pand	%xmm0, %xmm5
	paddd	%xmm5, %xmm3
	movups	%xmm3, 16(%rcx)
	movdqa	%xmm2, %xmm3
	psrld	$31, %xmm2
	pand	%xmm0, %xmm3
	paddd	%xmm3, %xmm2
	movdqa	%xmm1, %xmm3
	pand	%xmm0, %xmm3
	psrld	$31, %xmm1
	paddd	%xmm1, %xmm3
	movdqa	%xmm2, %xmm1
	pslld	$15, %xmm1
	psrld	$16, %xmm2
	pand	%xmm15, %xmm1
	por	%xmm2, %xmm1
	movdqa	%xmm3, %xmm2
	pslld	$15, %xmm2
	psrld	$16, %xmm3
	pand	%xmm15, %xmm2
	por	%xmm3, %xmm2
	movdqa	%xmm2, %xmm3
	paddd	%xmm1, %xmm3
	pxor	%xmm0, %xmm1
	paddd	%xmm2, %xmm1
	movdqa	%xmm3, %xmm4
	movdqa	%xmm1, %xmm2
	psrld	$31, %xmm3
	pand	%xmm0, %xmm4
	pand	%xmm0, %xmm2
	psrld	$31, %xmm1
	paddd	%xmm4, %xmm3
	paddd	%xmm2, %xmm1
	movups	%xmm3, 16(%rsi)
	movups	%xmm1, 16(%rax)
	cmpl	$2, %ebx
	je	.L166
	movdqu	32(%rcx), %xmm2
	movdqu	32(%rdx), %xmm6
	movdqu	32(%rax), %xmm1
	movdqu	32(%rsi), %xmm4
	movdqa	%xmm2, %xmm3
	pandn	%xmm0, %xmm2
	paddd	%xmm6, %xmm3
	paddd	%xmm6, %xmm2
	movdqa	%xmm3, %xmm5
	psrld	$31, %xmm3
	pand	%xmm0, %xmm5
	paddd	%xmm5, %xmm3
	movups	%xmm3, 32(%rdx)
	movdqa	%xmm1, %xmm3
	pandn	%xmm0, %xmm1
	paddd	%xmm4, %xmm3
	paddd	%xmm4, %xmm1
	movdqa	%xmm3, %xmm5
	pand	%xmm0, %xmm3
	psrld	$31, %xmm5
	paddd	%xmm5, %xmm3
	movups	%xmm3, 32(%rcx)
	movdqa	%xmm2, %xmm3
	pand	%xmm0, %xmm2
	psrld	$31, %xmm3
	paddd	%xmm3, %xmm2
	movdqa	%xmm1, %xmm3
	pand	%xmm0, %xmm1
	psrld	$31, %xmm3
	paddd	%xmm1, %xmm3
	movdqa	%xmm2, %xmm1
	pslld	$15, %xmm1
	psrld	$16, %xmm2
	pand	%xmm15, %xmm1
	por	%xmm2, %xmm1
	movdqa	%xmm3, %xmm2
	pslld	$15, %xmm2
	psrld	$16, %xmm3
	pand	%xmm15, %xmm2
	por	%xmm3, %xmm2
	movdqa	%xmm1, %xmm3
	pxor	%xmm0, %xmm1
	paddd	%xmm2, %xmm3
	paddd	%xmm2, %xmm1
	movdqa	%xmm3, %xmm4
	movdqa	%xmm1, %xmm2
	pand	%xmm0, %xmm3
	psrld	$31, %xmm4
	psrld	$31, %xmm2
	pand	%xmm0, %xmm1
	paddd	%xmm4, %xmm3
	paddd	%xmm2, %xmm1
	movups	%xmm3, 32(%rsi)
	movups	%xmm1, 32(%rax)
	cmpl	$4, %ebx
	jne	.L166
	movdqu	48(%rdx), %xmm6
	movdqu	48(%rcx), %xmm2
	movdqu	48(%rsi), %xmm4
	movdqu	48(%rax), %xmm1
	movdqa	%xmm6, %xmm3
	paddd	%xmm2, %xmm3
	pandn	%xmm0, %xmm2
	movdqa	%xmm3, %xmm5
	pand	%xmm0, %xmm3
	paddd	%xmm6, %xmm2
	psrld	$31, %xmm5
	paddd	%xmm5, %xmm3
	movups	%xmm3, 48(%rdx)
	movdqa	%xmm4, %xmm3
	paddd	%xmm1, %xmm3
	pandn	%xmm0, %xmm1
	movdqa	%xmm3, %xmm5
	pand	%xmm0, %xmm3
	paddd	%xmm4, %xmm1
	psrld	$31, %xmm5
	paddd	%xmm5, %xmm3
	movups	%xmm3, 48(%rcx)
	movdqa	%xmm2, %xmm3
	pand	%xmm0, %xmm2
	psrld	$31, %xmm3
	paddd	%xmm3, %xmm2
	movdqa	%xmm1, %xmm3
	pand	%xmm0, %xmm1
	psrld	$31, %xmm3
	paddd	%xmm1, %xmm3
	movdqa	%xmm2, %xmm1
	pslld	$15, %xmm1
	psrld	$16, %xmm2
	pand	%xmm15, %xmm1
	por	%xmm2, %xmm1
	movdqa	%xmm3, %xmm2
	pslld	$15, %xmm2
	psrld	$16, %xmm3
	pand	%xmm15, %xmm2
	por	%xmm3, %xmm2
	movdqa	%xmm1, %xmm3
	pxor	%xmm0, %xmm1
	paddd	%xmm2, %xmm3
	paddd	%xmm2, %xmm1
	movdqa	%xmm3, %xmm4
	movdqa	%xmm1, %xmm2
	pand	%xmm0, %xmm3
	psrld	$31, %xmm4
	psrld	$31, %xmm2
	pand	%xmm0, %xmm1
	paddd	%xmm4, %xmm3
	paddd	%xmm2, %xmm1
	movups	%xmm3, 48(%rsi)
	movups	%xmm1, 48(%rax)
	.p2align 4,,10
	.p2align 3
.L171:
	movl	-120(%rsp), %eax
	movl	-36(%rsp), %esi
	addl	%esi, %eax
	movl	%eax, -64(%rsp)
.L169:
	leaq	(%r12,%r12,2), %rbx
	leaq	0(%r13,%r15), %rax
	leaq	(%rbx,%r13), %rsi
	leaq	(%rax,%rbx), %rbp
	movq	%rsi, -80(%rsp)
	leaq	(%rbx,%r15), %rax
	salq	$2, %rsi
	leaq	16(,%rbp,4), %r9
	leaq	0(,%rax,4), %rcx
	leaq	16(%rsi), %r11
	movq	%rax, -56(%rsp)
	leaq	0(,%rbp,4), %rax
	leaq	0(,%rbx,4), %rdi
	cmpq	%rax, %r11
	leaq	16(%rdi), %rdx
	leaq	16(%rcx), %r8
	setle	%r10b
	cmpq	%r9, %rsi
	movq	%rdx, -120(%rsp)
	setge	%dl
	orl	%edx, %r10d
	cmpq	%r8, %rax
	setge	%dl
	cmpq	%r9, %rcx
	setge	-104(%rsp)
	orb	-104(%rsp), %dl
	andl	%r10d, %edx
	cmpq	%rax, -120(%rsp)
	setle	%r10b
	cmpq	%r9, %rdi
	setge	%r9b
	orl	%r10d, %r9d
	andl	%edx, %r9d
	cmpq	%rcx, %r11
	setle	%dl
	cmpq	%r8, %rsi
	setge	%r10b
	orl	%r10d, %edx
	andl	%r9d, %edx
	cmpq	%rcx, -120(%rsp)
	setle	%r9b
	cmpq	%r8, %rdi
	setge	%r8b
	orl	%r8d, %r9d
	testb	%r9b, %dl
	je	.L304
	cmpq	%rsi, -120(%rsp)
	setle	%dl
	cmpq	%r11, %rdi
	setge	%r8b
	orb	%dl, %r8b
	je	.L304
	movl	-72(%rsp), %r9d
	leal	-1(%r9), %edx
	movl	%r9d, %r8d
	cmpl	$2, %edx
	jbe	.L197
	addq	%r14, %rcx
	leaq	(%r14,%rdi), %rdx
	addq	%r14, %rax
	addq	%r14, %rsi
	movdqu	(%rcx), %xmm6
	movdqu	(%rdx), %xmm2
	shrl	$2, %r9d
	movdqu	(%rax), %xmm1
	movdqu	(%rsi), %xmm5
	movdqa	%xmm6, %xmm3
	paddd	%xmm2, %xmm3
	pandn	%xmm0, %xmm2
	movdqa	%xmm3, %xmm4
	psrld	$31, %xmm3
	paddd	%xmm6, %xmm2
	pand	%xmm0, %xmm4
	paddd	%xmm4, %xmm3
	movups	%xmm3, (%rdx)
	movdqa	%xmm1, %xmm3
	pandn	%xmm0, %xmm1
	paddd	%xmm5, %xmm3
	paddd	%xmm5, %xmm1
	movdqa	%xmm3, %xmm4
	psrld	$31, %xmm3
	pand	%xmm0, %xmm4
	paddd	%xmm4, %xmm3
	movdqa	%xmm2, %xmm4
	pand	%xmm0, %xmm2
	movups	%xmm3, (%rcx)
	psrld	$31, %xmm4
	movdqa	%xmm1, %xmm3
	pand	%xmm0, %xmm1
	paddd	%xmm2, %xmm4
	psrld	$31, %xmm3
	paddd	%xmm1, %xmm3
	movdqa	%xmm4, %xmm2
	pslld	$15, %xmm2
	psrld	$16, %xmm4
	movdqa	%xmm3, %xmm1
	pslld	$15, %xmm1
	psrld	$16, %xmm3
	pand	%xmm15, %xmm2
	por	%xmm4, %xmm2
	pand	%xmm15, %xmm1
	por	%xmm3, %xmm1
	movdqa	%xmm2, %xmm3
	paddd	%xmm1, %xmm3
	pxor	%xmm0, %xmm1
	paddd	%xmm2, %xmm1
	movdqa	%xmm3, %xmm4
	pand	%xmm0, %xmm3
	movdqa	%xmm1, %xmm2
	psrld	$31, %xmm4
	pand	%xmm0, %xmm2
	psrld	$31, %xmm1
	paddd	%xmm4, %xmm3
	paddd	%xmm2, %xmm1
	movups	%xmm3, (%rsi)
	movups	%xmm1, (%rax)
	cmpl	$1, %r9d
	je	.L176
	movdqu	16(%rdx), %xmm2
	movdqu	16(%rcx), %xmm6
	movdqu	16(%rsi), %xmm5
	movdqu	16(%rax), %xmm1
	movdqa	%xmm2, %xmm3
	pandn	%xmm0, %xmm2
	paddd	%xmm6, %xmm3
	paddd	%xmm6, %xmm2
	movdqa	%xmm3, %xmm4
	pand	%xmm0, %xmm3
	psrld	$31, %xmm4
	paddd	%xmm4, %xmm3
	movups	%xmm3, 16(%rdx)
	movdqa	%xmm5, %xmm3
	paddd	%xmm1, %xmm3
	pandn	%xmm0, %xmm1
	movdqa	%xmm3, %xmm4
	pand	%xmm0, %xmm3
	paddd	%xmm5, %xmm1
	psrld	$31, %xmm4
	paddd	%xmm4, %xmm3
	movdqa	%xmm2, %xmm4
	pand	%xmm0, %xmm2
	movups	%xmm3, 16(%rcx)
	psrld	$31, %xmm4
	movdqa	%xmm1, %xmm3
	pand	%xmm0, %xmm1
	paddd	%xmm2, %xmm4
	psrld	$31, %xmm3
	paddd	%xmm1, %xmm3
	movdqa	%xmm4, %xmm2
	pslld	$15, %xmm2
	psrld	$16, %xmm4
	movdqa	%xmm3, %xmm1
	pslld	$15, %xmm1
	psrld	$16, %xmm3
	pand	%xmm15, %xmm2
	por	%xmm4, %xmm2
	pand	%xmm15, %xmm1
	por	%xmm3, %xmm1
	movdqa	%xmm2, %xmm3
	paddd	%xmm1, %xmm3
	pxor	%xmm0, %xmm1
	paddd	%xmm2, %xmm1
	movdqa	%xmm3, %xmm4
	pand	%xmm0, %xmm3
	movdqa	%xmm1, %xmm2
	psrld	$31, %xmm4
	pand	%xmm0, %xmm1
	psrld	$31, %xmm2
	paddd	%xmm4, %xmm3
	paddd	%xmm2, %xmm1
	movups	%xmm3, 16(%rsi)
	movups	%xmm1, 16(%rax)
	cmpl	$2, %r9d
	je	.L176
	movdqu	32(%rdx), %xmm2
	movdqu	32(%rcx), %xmm6
	movdqu	32(%rsi), %xmm5
	movdqu	32(%rax), %xmm1
	movdqa	%xmm2, %xmm3
	pandn	%xmm0, %xmm2
	paddd	%xmm6, %xmm3
	paddd	%xmm6, %xmm2
	movdqa	%xmm3, %xmm4
	pand	%xmm0, %xmm3
	psrld	$31, %xmm4
	paddd	%xmm4, %xmm3
	movups	%xmm3, 32(%rdx)
	movdqa	%xmm5, %xmm3
	paddd	%xmm1, %xmm3
	pandn	%xmm0, %xmm1
	movdqa	%xmm3, %xmm4
	psrld	$31, %xmm3
	paddd	%xmm5, %xmm1
	pand	%xmm0, %xmm4
	paddd	%xmm4, %xmm3
	movdqa	%xmm2, %xmm4
	movups	%xmm3, 32(%rcx)
	psrld	$31, %xmm2
	movdqa	%xmm1, %xmm3
	pand	%xmm0, %xmm4
	paddd	%xmm2, %xmm4
	psrld	$31, %xmm3
	pand	%xmm0, %xmm1
	paddd	%xmm1, %xmm3
	movdqa	%xmm4, %xmm2
	pslld	$15, %xmm2
	psrld	$16, %xmm4
	movdqa	%xmm3, %xmm1
	pslld	$15, %xmm1
	psrld	$16, %xmm3
	pand	%xmm15, %xmm2
	por	%xmm4, %xmm2
	pand	%xmm15, %xmm1
	por	%xmm3, %xmm1
	movdqa	%xmm2, %xmm3
	paddd	%xmm1, %xmm3
	pxor	%xmm0, %xmm1
	paddd	%xmm2, %xmm1
	movdqa	%xmm3, %xmm4
	pand	%xmm0, %xmm3
	movdqa	%xmm1, %xmm2
	psrld	$31, %xmm4
	pand	%xmm0, %xmm2
	psrld	$31, %xmm1
	paddd	%xmm4, %xmm3
	paddd	%xmm2, %xmm1
	movups	%xmm3, 32(%rsi)
	movups	%xmm1, 32(%rax)
	cmpl	$4, %r9d
	jne	.L176
	movdqu	48(%rdx), %xmm2
	movdqu	48(%rcx), %xmm6
	movdqu	48(%rsi), %xmm5
	movdqu	48(%rax), %xmm1
	movdqa	%xmm2, %xmm3
	pandn	%xmm0, %xmm2
	paddd	%xmm6, %xmm3
	paddd	%xmm6, %xmm2
	movdqa	%xmm3, %xmm4
	pand	%xmm0, %xmm3
	psrld	$31, %xmm4
	paddd	%xmm4, %xmm3
	movups	%xmm3, 48(%rdx)
	movdqa	%xmm5, %xmm3
	paddd	%xmm1, %xmm3
	pandn	%xmm0, %xmm1
	movdqa	%xmm3, %xmm4
	pand	%xmm0, %xmm3
	paddd	%xmm5, %xmm1
	psrld	$31, %xmm4
	paddd	%xmm4, %xmm3
	movdqa	%xmm2, %xmm4
	pand	%xmm0, %xmm2
	movups	%xmm3, 48(%rcx)
	psrld	$31, %xmm4
	movdqa	%xmm1, %xmm3
	pand	%xmm0, %xmm1
	paddd	%xmm2, %xmm4
	psrld	$31, %xmm3
	paddd	%xmm1, %xmm3
	movdqa	%xmm4, %xmm2
	pslld	$15, %xmm2
	psrld	$16, %xmm4
	movdqa	%xmm3, %xmm1
	pslld	$15, %xmm1
	psrld	$16, %xmm3
	pand	%xmm15, %xmm2
	por	%xmm4, %xmm2
	pand	%xmm15, %xmm1
	por	%xmm3, %xmm1
	movdqa	%xmm2, %xmm3
	paddd	%xmm1, %xmm3
	pxor	%xmm0, %xmm1
	paddd	%xmm2, %xmm1
	movdqa	%xmm3, %xmm4
	pand	%xmm0, %xmm3
	movdqa	%xmm1, %xmm2
	psrld	$31, %xmm4
	pand	%xmm0, %xmm1
	psrld	$31, %xmm2
	paddd	%xmm4, %xmm3
	paddd	%xmm2, %xmm1
	movups	%xmm3, 48(%rsi)
	movups	%xmm1, 48(%rax)
	.p2align 4,,10
	.p2align 3
.L181:
	cmpl	$4, 40(%rsp)
	jle	.L143
	leaq	32+omegas128(%rip), %rax
	movq	%r14, -48(%rsp)
	movq	%rax, -104(%rsp)
	leaq	0(,%r12,4), %rax
	movq	%rax, -8(%rsp)
	movslq	-64(%rsp), %rax
	leaq	(%rax,%r12), %rdx
	leaq	(%r14,%rdx,4), %rdi
	leaq	0(%r13,%r12), %rdx
	addq	%rax, %rdx
	leaq	0(,%rdx,4), %rsi
	leaq	4(%r12,%r13), %rdx
	addq	%rax, %rdx
	leaq	4(%rax,%r12), %rax
	leaq	(%rsi,%r14), %rcx
	movq	%rsi, -80(%rsp)
	leaq	0(,%rax,4), %r12
	movl	40(%rsp), %eax
	leaq	0(,%rdx,4), %rsi
	leaq	40+omegas128(%rip), %rdx
	movq	%rsi, -120(%rsp)
	subl	$5, %eax
	leaq	(%rdx,%rax,8), %rax
	movq	%rax, -16(%rsp)
	movl	-72(%rsp), %eax
	leal	-1(%rax), %esi
	movl	%esi, (%rsp)
	leaq	0(,%r15,4), %rsi
	movq	%rsi, -32(%rsp)
	subq	%r14, %rsi
	movq	%rdi, %r14
	movq	%rsi, -24(%rsp)
	movl	%eax, %esi
	shrl	$2, %eax
	salq	$4, %rax
	movq	%rax, 16(%rsp)
	movl	%esi, %eax
	andl	$-4, %esi
	leal	1(%rsi), %ebx
	andl	$3, %eax
	movl	%esi, 32(%rsp)
	addl	$2, %esi
	movl	%ebx, 28(%rsp)
	movl	%eax, 24(%rsp)
	movl	%esi, 36(%rsp)
	movslq	-68(%rsp), %rsi
	.p2align 4,,10
	.p2align 3
.L189:
	movq	-104(%rsp), %rax
	movl	-36(%rsp), %edx
	addl	%edx, -64(%rsp)
	cmpl	$2, (%rsp)
	movl	4(%rax), %edx
	movd	(%rax), %xmm8
	movq	%rdx, %xmm7
	movd	%xmm8, %eax
	jbe	.L198
	movq	-24(%rsp), %rdi
	movq	-120(%rsp), %rbx
	leaq	-16(%r12), %rbp
	movq	-32(%rsp), %r8
	leaq	(%rdi,%rcx), %r11
	leaq	(%rdi,%r14), %r9
	leaq	(%rbx,%r8), %r10
	addq	%r12, %r8
	cmpq	%r11, %rbx
	setle	%bl
	cmpq	%r10, -80(%rsp)
	setge	%dil
	orl	%edi, %ebx
	cmpq	%r11, %r8
	setle	%dil
	cmpq	%r10, %r9
	setge	-56(%rsp)
	orb	-56(%rsp), %dil
	andl	%ebx, %edi
	cmpq	%r11, %r12
	movq	-120(%rsp), %rbx
	setle	%r11b
	cmpq	%r10, %rbp
	setge	%r10b
	orl	%r10d, %r11d
	andl	%r11d, %edi
	movq	-80(%rsp), %r11
	cmpq	%r9, %rbx
	setle	%r10b
	cmpq	%r11, %r8
	setle	%r11b
	orl	%r11d, %r10d
	andl	%edi, %r10d
	cmpq	%r12, %r9
	setge	%dil
	cmpq	%rbp, %r8
	setle	%r8b
	orl	%r8d, %edi
	testb	%dil, %r10b
	je	.L198
	cmpq	%r12, -80(%rsp)
	setge	%dil
	cmpq	%rbp, %rbx
	setle	%r8b
	orb	%dil, %r8b
	je	.L198
	pshufd	$0, %xmm8, %xmm8
	pshufd	$0, %xmm7, %xmm7
	movq	-32(%rsp), %rbx
	xorl	%edi, %edi
	movdqa	%xmm8, %xmm12
	movdqa	%xmm7, %xmm11
	movq	16(%rsp), %r10
	punpckldq	%xmm8, %xmm12
	punpckldq	%xmm7, %xmm11
	leaq	(%rbx,%r14), %r9
	leaq	(%rbx,%rcx), %r8
	punpckhdq	%xmm8, %xmm8
	punpckhdq	%xmm7, %xmm7
.L190:
	movdqu	(%r9,%rdi), %xmm1
	movdqu	(%r14,%rdi), %xmm5
	movdqa	%xmm7, %xmm10
	movdqu	(%r8,%rdi), %xmm3
	movdqu	(%rcx,%rdi), %xmm4
	movdqa	%xmm1, %xmm2
	pandn	%xmm0, %xmm1
	paddd	%xmm5, %xmm2
	paddd	%xmm5, %xmm1
	movdqa	%xmm2, %xmm6
	psrld	$31, %xmm2
	pand	%xmm0, %xmm6
	paddd	%xmm6, %xmm2
	movups	%xmm2, (%r14,%rdi)
	movdqa	%xmm3, %xmm2
	pandn	%xmm0, %xmm3
	paddd	%xmm4, %xmm2
	paddd	%xmm4, %xmm3
	movdqa	%xmm8, %xmm4
	movdqa	%xmm2, %xmm6
	psrld	$31, %xmm2
	pand	%xmm0, %xmm6
	paddd	%xmm6, %xmm2
	movups	%xmm2, (%r9,%rdi)
	movdqa	%xmm1, %xmm2
	psrld	$31, %xmm1
	pand	%xmm0, %xmm2
	paddd	%xmm1, %xmm2
	movdqa	%xmm3, %xmm1
	movdqa	%xmm2, %xmm5
	psrld	$31, %xmm3
	pand	%xmm0, %xmm1
	paddd	%xmm3, %xmm1
	punpckldq	%xmm2, %xmm5
	movdqa	%xmm12, %xmm3
	punpckhdq	%xmm2, %xmm2
	pmuludq	%xmm5, %xmm3
	pmuludq	%xmm2, %xmm4
	pmuludq	%xmm11, %xmm5
	pmuludq	%xmm7, %xmm2
	movdqa	%xmm3, %xmm6
	movdqa	%xmm4, %xmm9
	shufps	$136, %xmm4, %xmm3
	psrlq	$31, %xmm6
	pand	%xmm0, %xmm3
	psrlq	$31, %xmm9
	movdqa	%xmm1, %xmm4
	punpckldq	%xmm1, %xmm4
	shufps	$136, %xmm9, %xmm6
	punpckhdq	%xmm1, %xmm1
	paddd	%xmm3, %xmm6
	movdqa	%xmm11, %xmm3
	pmuludq	%xmm1, %xmm10
	pmuludq	%xmm4, %xmm3
	pmuludq	%xmm8, %xmm1
	pmuludq	%xmm12, %xmm4
	movdqa	%xmm10, %xmm13
	movdqa	%xmm3, %xmm9
	psrlq	$31, %xmm13
	shufps	$136, %xmm10, %xmm3
	pand	%xmm0, %xmm3
	psrlq	$31, %xmm9
	shufps	$136, %xmm13, %xmm9
	paddd	%xmm9, %xmm3
	movdqa	%xmm6, %xmm9
	pand	%xmm0, %xmm9
	psrld	$31, %xmm6
	paddd	%xmm9, %xmm6
	movdqa	%xmm3, %xmm9
	pand	%xmm0, %xmm9
	psrld	$31, %xmm3
	paddd	%xmm9, %xmm3
	paddd	%xmm6, %xmm3
	movdqa	%xmm3, %xmm6
	psrld	$31, %xmm3
	pand	%xmm0, %xmm6
	paddd	%xmm6, %xmm3
	movdqa	%xmm2, %xmm6
	movups	%xmm3, (%rcx,%rdi)
	movdqa	%xmm5, %xmm3
	psrlq	$31, %xmm6
	shufps	$136, %xmm2, %xmm5
	psrlq	$31, %xmm3
	pand	%xmm0, %xmm5
	movdqa	%xmm4, %xmm2
	shufps	$136, %xmm1, %xmm4
	shufps	$136, %xmm6, %xmm3
	paddd	%xmm5, %xmm3
	movdqa	%xmm1, %xmm5
	pand	%xmm0, %xmm4
	psrlq	$31, %xmm2
	psrlq	$31, %xmm5
	movdqa	%xmm3, %xmm1
	pand	%xmm0, %xmm1
	psrld	$31, %xmm3
	shufps	$136, %xmm5, %xmm2
	paddd	%xmm4, %xmm2
	paddd	%xmm3, %xmm1
	movdqa	%xmm2, %xmm3
	pand	%xmm0, %xmm3
	psrld	$31, %xmm2
	pandn	%xmm0, %xmm1
	paddd	%xmm3, %xmm2
	paddd	%xmm2, %xmm1
	movdqa	%xmm1, %xmm2
	psrld	$31, %xmm1
	pand	%xmm0, %xmm2
	paddd	%xmm2, %xmm1
	movups	%xmm1, (%r8,%rdi)
	addq	$16, %rdi
	cmpq	%rdi, %r10
	jne	.L190
	movl	24(%rsp), %edi
	testl	%edi, %edi
	je	.L192
	movl	32(%rsp), %ebx
	movl	-64(%rsp), %edi
	movq	-48(%rsp), %r9
	movl	-68(%rsp), %r10d
	addl	%ebx, %edi
	movslq	%edi, %r8
	leaq	(%r9,%r8,4), %rbp
	movl	0(%rbp), %ebx
	movl	%ebx, 4(%rsp)
	movl	-72(%rsp), %ebx
	leal	(%rbx,%rdi), %r8d
	movslq	%r8d, %r8
	leaq	(%r9,%r8,4), %r11
	movl	(%r11), %r8d
	movl	%r8d, -56(%rsp)
	movl	-40(%rsp), %r8d
	addl	%edi, %r8d
	addl	%r10d, %edi
	movl	4(%rsp), %r10d
	movslq	%r8d, %r8
	movslq	%edi, %rdi
	leaq	(%r9,%r8,4), %rbx
	leaq	(%r9,%rdi,4), %rdi
	movl	(%rbx), %r8d
	movq	%rdi, 8(%rsp)
	movl	(%rdi), %edi
	leal	(%r10,%r8), %r9d
	notl	%r8d
	movl	%r9d, %r10d
	andl	$2147483647, %r9d
	andl	$2147483647, %r8d
	shrl	$31, %r10d
	addl	%r9d, %r10d
	movl	%r10d, 0(%rbp)
	movl	-56(%rsp), %ebp
	leal	0(%rbp,%rdi), %r9d
	movl	4(%rsp), %ebp
	notl	%edi
	andl	$2147483647, %edi
	movl	%r9d, %r10d
	andl	$2147483647, %r9d
	addl	%ebp, %r8d
	movl	-56(%rsp), %ebp
	shrl	$31, %r10d
	addl	%r9d, %r10d
	addl	%ebp, %edi
	movl	%r8d, %ebp
	shrl	$31, %r8d
	movl	%r10d, (%rbx)
	andl	$2147483647, %ebp
	addl	%ebp, %r8d
	movq	%rax, %rbp
	movl	%r8d, %r8d
	imulq	%r8, %rbp
	movq	%rbp, %r10
	shrq	$31, %rbp
	andl	$2147483647, %r10d
	addl	%ebp, %r10d
	movl	%edi, %ebp
	andl	$2147483647, %edi
	shrl	$31, %ebp
	addl	%ebp, %edi
	movq	%rdx, %rbp
	movl	%edi, %edi
	imulq	%rdi, %rbp
	movq	%rbp, %r9
	shrq	$31, %rbp
	andl	$2147483647, %r9d
	addl	%ebp, %r9d
	movl	%r10d, %ebp
	shrl	$31, %r10d
	movl	%r9d, %ebx
	andl	$2147483647, %ebp
	andl	$2147483647, %ebx
	addl	%ebp, %ebx
	addl	%ebx, %r10d
	shrl	$31, %r9d
	imulq	%rdx, %r8
	movl	-72(%rsp), %ebx
	addl	%r10d, %r9d
	imulq	%rax, %rdi
	movl	%r9d, %ebp
	andl	$2147483647, %r9d
	shrl	$31, %ebp
	movl	%ebp, %r10d
	movq	%r8, %rbp
	andl	$2147483647, %r8d
	shrq	$31, %rbp
	addl	%r9d, %r10d
	addl	%ebp, %r8d
	movq	%rdi, %rbp
	andl	$2147483647, %edi
	movl	%r10d, (%r11)
	shrq	$31, %rbp
	movq	8(%rsp), %r11
	addl	%ebp, %edi
	movl	%r8d, %ebp
	andl	$2147483647, %r8d
	shrl	$31, %ebp
	addl	%ebp, %r8d
	movl	%edi, %ebp
	andl	$2147483647, %edi
	shrl	$31, %ebp
	notl	%r8d
	andl	$2147483647, %r8d
	addl	%ebp, %edi
	addl	%r8d, %edi
	movl	%edi, %ebp
	shrl	$31, %edi
	andl	$2147483647, %ebp
	movl	%ebp, %r8d
	addl	%edi, %r8d
	movl	%r8d, (%r11)
	movl	28(%rsp), %r11d
	cmpl	%r11d, %ebx
	jle	.L192
	movl	-64(%rsp), %edi
	movq	-48(%rsp), %r9
	movl	-68(%rsp), %r10d
	addl	%r11d, %edi
	movslq	%edi, %r8
	leaq	(%r9,%r8,4), %rbp
	leal	(%rbx,%rdi), %r8d
	movl	0(%rbp), %r11d
	movslq	%r8d, %r8
	movl	%r11d, -56(%rsp)
	leaq	(%r9,%r8,4), %r11
	movl	-40(%rsp), %r8d
	movl	(%r11), %ebx
	addl	%edi, %r8d
	addl	%r10d, %edi
	movl	-56(%rsp), %r10d
	movslq	%r8d, %r8
	movl	%ebx, 4(%rsp)
	movslq	%edi, %rdi
	leaq	(%r9,%r8,4), %rbx
	leaq	(%r9,%rdi,4), %rdi
	movl	(%rbx), %r8d
	movq	%rdi, 8(%rsp)
	movl	(%rdi), %edi
	leal	(%r10,%r8), %r9d
	notl	%r8d
	movl	%r9d, %r10d
	andl	$2147483647, %r9d
	andl	$2147483647, %r8d
	shrl	$31, %r10d
	addl	%r9d, %r10d
	movl	%r10d, 0(%rbp)
	movl	4(%rsp), %ebp
	leal	0(%rbp,%rdi), %r9d
	notl	%edi
	movl	%r9d, %ebp
	andl	$2147483647, %edi
	andl	$2147483647, %r9d
	shrl	$31, %ebp
	movl	%ebp, %r10d
	movl	-56(%rsp), %ebp
	addl	%r9d, %r10d
	addl	%ebp, %r8d
	movl	4(%rsp), %ebp
	movl	%r10d, (%rbx)
	addl	%ebp, %edi
	movl	%r8d, %ebp
	andl	$2147483647, %r8d
	shrl	$31, %ebp
	leal	0(%rbp,%r8), %r10d
	movq	%rax, %rbp
	imulq	%r10, %rbp
	movq	%rbp, %r9
	shrq	$31, %rbp
	andl	$2147483647, %r9d
	addl	%ebp, %r9d
	movl	%edi, %ebp
	andl	$2147483647, %edi
	shrl	$31, %ebp
	movl	%r9d, %ebx
	shrl	$31, %r9d
	addl	%ebp, %edi
	movq	%rdx, %rbp
	andl	$2147483647, %ebx
	movl	%edi, %edi
	imulq	%rdi, %rbp
	movq	%rbp, %r8
	shrq	$31, %rbp
	andl	$2147483647, %r8d
	addl	%ebp, %r8d
	movl	%r8d, %ebp
	shrl	$31, %r8d
	andl	$2147483647, %ebp
	addl	%ebp, %ebx
	addl	%ebx, %r9d
	movl	36(%rsp), %ebx
	addl	%r9d, %r8d
	imulq	%rdx, %r10
	movl	%r8d, %ebp
	imulq	%rax, %rdi
	shrl	$31, %r8d
	andl	$2147483647, %ebp
	movl	%ebp, %r9d
	movl	%r10d, %ebp
	shrq	$31, %r10
	andl	$2147483647, %ebp
	addl	%r8d, %r9d
	movl	%ebp, %r8d
	movl	%edi, %ebp
	shrq	$31, %rdi
	movl	%r9d, (%r11)
	addl	%r10d, %r8d
	andl	$2147483647, %ebp
	movq	8(%rsp), %r11
	movl	-72(%rsp), %r10d
	addl	%ebp, %edi
	movl	%r8d, %ebp
	shrl	$31, %r8d
	andl	$2147483647, %ebp
	addl	%ebp, %r8d
	movl	%edi, %ebp
	shrl	$31, %edi
	andl	$2147483647, %ebp
	notl	%r8d
	andl	$2147483647, %r8d
	addl	%ebp, %edi
	addl	%r8d, %edi
	movl	%edi, %ebp
	shrl	$31, %edi
	andl	$2147483647, %ebp
	movl	%ebp, %r8d
	addl	%edi, %r8d
	movl	%r8d, (%r11)
	cmpl	%ebx, %r10d
	jle	.L192
	movl	-64(%rsp), %edi
	movq	-48(%rsp), %r9
	addl	%ebx, %edi
	movslq	%edi, %r8
	leaq	(%r9,%r8,4), %r11
	leal	(%r10,%rdi), %r8d
	movl	-68(%rsp), %r10d
	movslq	%r8d, %r8
	movl	(%r11), %ebx
	leaq	(%r9,%r8,4), %rbp
	movl	0(%rbp), %r8d
	movq	%rbp, 8(%rsp)
	movl	%r8d, 4(%rsp)
	movl	-40(%rsp), %r8d
	addl	%edi, %r8d
	addl	%r10d, %edi
	movslq	%r8d, %r8
	movslq	%edi, %rdi
	leaq	(%r9,%r8,4), %rbp
	leaq	(%r9,%rdi,4), %rdi
	movl	0(%rbp), %r8d
	movq	%rdi, -56(%rsp)
	movl	(%rdi), %edi
	leal	(%r8,%rbx), %r9d
	notl	%r8d
	movl	%r9d, %r10d
	shrl	$31, %r9d
	andl	$2147483647, %r8d
	andl	$2147483647, %r10d
	addl	%ebx, %r8d
	addl	%r9d, %r10d
	movl	%r10d, (%r11)
	movl	4(%rsp), %r11d
	leal	(%rdi,%r11), %r9d
	notl	%edi
	movl	%r9d, %r10d
	shrl	$31, %r9d
	andl	$2147483647, %edi
	andl	$2147483647, %r10d
	addl	%r11d, %edi
	addl	%r9d, %r10d
	movl	%r8d, %r9d
	shrl	$31, %r8d
	andl	$2147483647, %r9d
	movl	%r10d, 0(%rbp)
	movq	8(%rsp), %rbp
	addl	%r9d, %r8d
	movl	%r8d, %r8d
	movq	%r8, %r9
	imulq	%rax, %r9
	movl	%r9d, %r10d
	shrq	$31, %r9
	andl	$2147483647, %r10d
	addl	%r9d, %r10d
	movl	%edi, %r9d
	shrl	$31, %edi
	andl	$2147483647, %r9d
	addl	%r9d, %edi
	movl	%edi, %edi
	movq	%rdi, %r11
	imulq	%rdx, %r11
	movl	%r11d, %r9d
	shrq	$31, %r11
	andl	$2147483647, %r9d
	addl	%r11d, %r9d
	movl	%r10d, %r11d
	shrl	$31, %r10d
	movl	%r9d, %ebx
	andl	$2147483647, %r11d
	andl	$2147483647, %ebx
	addl	%ebx, %r11d
	addl	%r11d, %r10d
	shrl	$31, %r9d
	imulq	%rdx, %r8
	imulq	%rax, %rdi
	addl	%r10d, %r9d
	movl	%r9d, %r10d
	shrl	$31, %r9d
	movl	%r8d, %edx
	shrq	$31, %r8
	andl	$2147483647, %r10d
	andl	$2147483647, %edx
	addl	%r9d, %r10d
	addl	%edx, %r8d
	movl	%edi, %edx
	shrq	$31, %rdi
	movl	%r10d, 0(%rbp)
	andl	$2147483647, %edx
	leal	(%rdx,%rdi), %eax
	movl	%r8d, %edx
	shrl	$31, %r8d
	movq	-56(%rsp), %rdi
	andl	$2147483647, %edx
	movl	%eax, %ebx
	shrl	$31, %eax
	addl	%r8d, %edx
	andl	$2147483647, %ebx
	notl	%edx
	addl	%ebx, %eax
	andl	$2147483647, %edx
	addl	%edx, %eax
	movl	%eax, %edx
	shrl	$31, %eax
	andl	$2147483647, %edx
	addl	%eax, %edx
	movl	%edx, (%rdi)
.L192:
	movq	-8(%rsp), %rax
	addq	$8, -104(%rsp)
	addq	%rax, -80(%rsp)
	movq	-104(%rsp), %rdx
	addq	%rax, -120(%rsp)
	addq	%rax, %r14
	addq	%rax, %rcx
	addq	%rax, %r12
	cmpq	%rdx, -16(%rsp)
	jne	.L189
	movq	-48(%rsp), %r14
.L143:
	sarl	40(%rsp)
	subl	$1, 44(%rsp)
	jne	.L183
	movq	48(%rsp), %rbx
	movdqu	256(%rbx), %xmm1
	movdqu	(%rbx), %xmm4
	movdqu	384(%rbx), %xmm7
	movdqu	128(%rbx), %xmm3
	movdqa	%xmm1, %xmm2
	pandn	%xmm0, %xmm1
	paddd	%xmm4, %xmm2
	paddd	%xmm4, %xmm1
	movdqu	400(%rbx), %xmm4
	movdqa	%xmm2, %xmm6
	pand	%xmm0, %xmm2
	psrld	$31, %xmm6
	movaps	%xmm4, -120(%rsp)
	paddd	%xmm2, %xmm6
	movdqa	%xmm7, %xmm2
	paddd	%xmm3, %xmm2
	movups	%xmm6, (%rbx)
	movd	%xmm6, %ecx
	movdqa	%xmm2, %xmm5
	pand	%xmm0, %xmm2
	psrld	$31, %xmm5
	paddd	%xmm2, %xmm5
	movdqa	%xmm1, %xmm2
	pand	%xmm0, %xmm1
	psrld	$31, %xmm2
	movups	%xmm5, 256(%rbx)
	movd	%xmm5, %edi
	paddd	%xmm2, %xmm1
	movups	%xmm1, 128(%rbx)
	movdqa	%xmm7, %xmm1
	movdqu	16(%rbx), %xmm7
	pandn	%xmm0, %xmm1
	paddd	%xmm3, %xmm1
	movdqu	144(%rbx), %xmm3
	movdqa	%xmm1, %xmm2
	psrld	$31, %xmm1
	pand	%xmm0, %xmm2
	paddd	%xmm2, %xmm1
	movdqa	%xmm7, %xmm2
	movups	%xmm1, 384(%rbx)
	movdqu	272(%rbx), %xmm1
	paddd	%xmm1, %xmm2
	pandn	%xmm0, %xmm1
	movdqa	%xmm2, %xmm4
	psrld	$31, %xmm2
	paddd	%xmm7, %xmm1
	movdqu	32(%rbx), %xmm7
	pand	%xmm0, %xmm4
	paddd	%xmm4, %xmm2
	movups	%xmm2, 16(%rbx)
	movdqa	-120(%rsp), %xmm2
	paddd	%xmm3, %xmm2
	movdqa	%xmm2, %xmm4
	pand	%xmm0, %xmm2
	psrld	$31, %xmm4
	paddd	%xmm4, %xmm2
	movdqu	416(%rbx), %xmm4
	movups	%xmm2, 272(%rbx)
	movdqa	%xmm1, %xmm2
	pand	%xmm0, %xmm1
	psrld	$31, %xmm2
	paddd	%xmm2, %xmm1
	movups	%xmm1, 144(%rbx)
	movdqa	-120(%rsp), %xmm1
	movaps	%xmm4, -120(%rsp)
	pandn	%xmm0, %xmm1
	paddd	%xmm3, %xmm1
	movdqu	160(%rbx), %xmm3
	movdqa	%xmm1, %xmm2
	pand	%xmm0, %xmm1
	psrld	$31, %xmm2
	paddd	%xmm2, %xmm1
	movdqa	%xmm7, %xmm2
	movups	%xmm1, 400(%rbx)
	movdqu	288(%rbx), %xmm1
	paddd	%xmm1, %xmm2
	pandn	%xmm0, %xmm1
	movdqa	%xmm2, %xmm4
	pand	%xmm0, %xmm2
	paddd	%xmm7, %xmm1
	movdqu	48(%rbx), %xmm7
	psrld	$31, %xmm4
	paddd	%xmm4, %xmm2
	movups	%xmm2, 32(%rbx)
	movdqa	-120(%rsp), %xmm2
	paddd	%xmm3, %xmm2
	movdqa	%xmm2, %xmm4
	pand	%xmm0, %xmm2
	psrld	$31, %xmm4
	paddd	%xmm4, %xmm2
	movups	%xmm2, 288(%rbx)
	movdqa	%xmm1, %xmm2
	pand	%xmm0, %xmm1
	psrld	$31, %xmm2
	paddd	%xmm2, %xmm1
	movups	%xmm1, 160(%rbx)
	movdqa	-120(%rsp), %xmm1
	pandn	%xmm0, %xmm1
	paddd	%xmm3, %xmm1
	movdqa	%xmm1, %xmm2
	pand	%xmm0, %xmm1
	psrld	$31, %xmm2
	paddd	%xmm2, %xmm1
	movdqa	%xmm7, %xmm2
	movups	%xmm1, 416(%rbx)
	movdqu	176(%rbx), %xmm3
	movdqu	304(%rbx), %xmm1
	movdqu	432(%rbx), %xmm4
	paddd	%xmm1, %xmm2
	pandn	%xmm0, %xmm1
	movaps	%xmm4, -120(%rsp)
	movdqa	%xmm2, %xmm4
	pand	%xmm0, %xmm2
	paddd	%xmm7, %xmm1
	psrld	$31, %xmm4
	movdqu	64(%rbx), %xmm7
	paddd	%xmm4, %xmm2
	movups	%xmm2, 48(%rbx)
	movdqa	-120(%rsp), %xmm2
	paddd	%xmm3, %xmm2
	movdqa	%xmm2, %xmm4
	pand	%xmm0, %xmm2
	psrld	$31, %xmm4
	paddd	%xmm4, %xmm2
	movdqu	448(%rbx), %xmm4
	movups	%xmm2, 304(%rbx)
	movdqa	%xmm1, %xmm2
	psrld	$31, %xmm1
	pand	%xmm0, %xmm2
	paddd	%xmm2, %xmm1
	movups	%xmm1, 176(%rbx)
	movdqa	-120(%rsp), %xmm1
	movaps	%xmm4, -120(%rsp)
	pandn	%xmm0, %xmm1
	paddd	%xmm3, %xmm1
	movdqu	192(%rbx), %xmm3
	movdqa	%xmm1, %xmm2
	pand	%xmm0, %xmm1
	psrld	$31, %xmm2
	paddd	%xmm2, %xmm1
	movdqa	%xmm7, %xmm2
	movups	%xmm1, 432(%rbx)
	movdqu	320(%rbx), %xmm1
	paddd	%xmm1, %xmm2
	pandn	%xmm0, %xmm1
	movdqa	%xmm2, %xmm4
	pand	%xmm0, %xmm2
	paddd	%xmm7, %xmm1
	movdqu	80(%rbx), %xmm7
	psrld	$31, %xmm4
	paddd	%xmm4, %xmm2
	movups	%xmm2, 64(%rbx)
	movdqa	-120(%rsp), %xmm2
	paddd	%xmm3, %xmm2
	movdqa	%xmm2, %xmm4
	pand	%xmm0, %xmm2
	psrld	$31, %xmm4
	paddd	%xmm4, %xmm2
	movdqu	464(%rbx), %xmm4
	movups	%xmm2, 320(%rbx)
	movdqa	%xmm1, %xmm2
	pand	%xmm0, %xmm1
	psrld	$31, %xmm2
	paddd	%xmm2, %xmm1
	movups	%xmm1, 192(%rbx)
	movdqa	-120(%rsp), %xmm1
	movaps	%xmm4, -120(%rsp)
	pandn	%xmm0, %xmm1
	paddd	%xmm3, %xmm1
	movdqu	208(%rbx), %xmm3
	movdqa	%xmm1, %xmm2
	pand	%xmm0, %xmm1
	psrld	$31, %xmm2
	paddd	%xmm2, %xmm1
	movdqa	%xmm7, %xmm2
	movups	%xmm1, 448(%rbx)
	movdqu	336(%rbx), %xmm1
	paddd	%xmm1, %xmm2
	pandn	%xmm0, %xmm1
	movdqa	%xmm2, %xmm4
	pand	%xmm0, %xmm2
	paddd	%xmm7, %xmm1
	psrld	$31, %xmm4
	paddd	%xmm4, %xmm2
	movups	%xmm2, 80(%rbx)
	movdqa	-120(%rsp), %xmm2
	paddd	%xmm3, %xmm2
	movdqa	%xmm2, %xmm4
	pand	%xmm0, %xmm2
	psrld	$31, %xmm4
	paddd	%xmm4, %xmm2
	movups	%xmm2, 336(%rbx)
	movdqa	%xmm1, %xmm2
	psrld	$31, %xmm1
	pand	%xmm0, %xmm2
	paddd	%xmm2, %xmm1
	movups	%xmm1, 208(%rbx)
	movdqa	-120(%rsp), %xmm1
	pandn	%xmm0, %xmm1
	paddd	%xmm3, %xmm1
	movdqa	%xmm1, %xmm2
	pand	%xmm0, %xmm1
	psrld	$31, %xmm2
	paddd	%xmm2, %xmm1
	movups	%xmm1, 464(%rbx)
	movdqu	96(%rbx), %xmm7
	movdqu	352(%rbx), %xmm1
	movdqu	480(%rbx), %xmm4
	movdqu	224(%rbx), %xmm3
	movdqa	%xmm7, %xmm2
	paddd	%xmm1, %xmm2
	movaps	%xmm4, -120(%rsp)
	pandn	%xmm0, %xmm1
	movdqa	%xmm2, %xmm4
	psrld	$31, %xmm2
	paddd	%xmm7, %xmm1
	movdqu	112(%rbx), %xmm7
	pand	%xmm0, %xmm4
	paddd	%xmm4, %xmm2
	movups	%xmm2, 96(%rbx)
	movdqa	-120(%rsp), %xmm2
	paddd	%xmm3, %xmm2
	movdqa	%xmm2, %xmm4
	pand	%xmm0, %xmm2
	psrld	$31, %xmm4
	paddd	%xmm4, %xmm2
	movdqu	496(%rbx), %xmm4
	movups	%xmm2, 352(%rbx)
	movdqa	%xmm1, %xmm2
	psrld	$31, %xmm1
	pand	%xmm0, %xmm2
	paddd	%xmm2, %xmm1
	movups	%xmm1, 224(%rbx)
	movdqa	-120(%rsp), %xmm1
	movaps	%xmm4, -120(%rsp)
	pandn	%xmm0, %xmm1
	paddd	%xmm3, %xmm1
	movdqu	240(%rbx), %xmm3
	movdqa	%xmm1, %xmm2
	pand	%xmm0, %xmm1
	psrld	$31, %xmm2
	paddd	%xmm2, %xmm1
	movups	%xmm1, 480(%rbx)
	movdqu	368(%rbx), %xmm1
	movdqa	%xmm1, %xmm2
	pandn	%xmm0, %xmm1
	paddd	%xmm7, %xmm2
	paddd	%xmm7, %xmm1
	movdqu	640(%rbx), %xmm7
	movdqa	%xmm2, %xmm4
	psrld	$31, %xmm2
	pand	%xmm0, %xmm4
	paddd	%xmm4, %xmm2
	movups	%xmm2, 112(%rbx)
	movdqa	-120(%rsp), %xmm2
	paddd	%xmm3, %xmm2
	movdqa	%xmm2, %xmm4
	psrld	$31, %xmm2
	pand	%xmm0, %xmm4
	paddd	%xmm4, %xmm2
	movdqu	896(%rbx), %xmm4
	movups	%xmm2, 368(%rbx)
	movdqa	%xmm1, %xmm2
	psrld	$31, %xmm1
	pand	%xmm0, %xmm2
	movaps	%xmm4, -104(%rsp)
	paddd	%xmm2, %xmm1
	movups	%xmm1, 240(%rbx)
	movdqa	-120(%rsp), %xmm1
	movaps	%xmm7, -120(%rsp)
	movdqu	768(%rbx), %xmm7
	movdqa	-120(%rsp), %xmm4
	pandn	%xmm0, %xmm1
	paddd	-104(%rsp), %xmm4
	paddd	%xmm3, %xmm1
	movdqu	512(%rbx), %xmm3
	movdqa	%xmm1, %xmm2
	psrld	$31, %xmm1
	pand	%xmm0, %xmm2
	paddd	%xmm2, %xmm1
	movups	%xmm1, 496(%rbx)
	movdqa	%xmm7, %xmm1
	paddd	%xmm3, %xmm1
	pandn	%xmm0, %xmm3
	movdqa	%xmm1, %xmm2
	psrld	$31, %xmm1
	paddd	%xmm7, %xmm3
	pand	%xmm0, %xmm2
	paddd	%xmm1, %xmm2
	movdqa	%xmm4, %xmm1
	pand	%xmm0, %xmm4
	psrld	$31, %xmm1
	movups	%xmm2, 512(%rbx)
	movd	%xmm2, %eax
	paddd	%xmm4, %xmm1
	leal	(%rax,%rcx), %esi
	notl	%eax
	movups	%xmm1, 768(%rbx)
	movdqa	-104(%rsp), %xmm4
	movl	%esi, %r8d
	movdqu	656(%rbx), %xmm7
	shrl	$31, %r8d
	andl	$2147483647, %esi
	movd	%xmm1, %edx
	andl	$2147483647, %eax
	pandn	%xmm0, %xmm4
	paddd	-120(%rsp), %xmm4
	addl	%r8d, %esi
	addl	%ecx, %eax
	movaps	%xmm7, -120(%rsp)
	movdqu	784(%rbx), %xmm7
	movdqa	%xmm4, %xmm8
	pand	%xmm0, %xmm4
	psrld	$31, %xmm8
	paddd	%xmm8, %xmm4
	movups	%xmm4, 640(%rbx)
	movdqa	%xmm3, %xmm4
	psrld	$31, %xmm3
	pand	%xmm0, %xmm4
	paddd	%xmm4, %xmm3
	movdqu	912(%rbx), %xmm4
	movups	%xmm3, 896(%rbx)
	movdqu	528(%rbx), %xmm3
	movaps	%xmm4, -104(%rsp)
	movdqa	%xmm7, %xmm4
	paddd	%xmm3, %xmm4
	pandn	%xmm0, %xmm3
	movdqa	%xmm4, %xmm8
	psrld	$31, %xmm4
	paddd	%xmm7, %xmm3
	movdqu	672(%rbx), %xmm7
	pand	%xmm0, %xmm8
	paddd	%xmm8, %xmm4
	movups	%xmm4, 528(%rbx)
	movdqa	-120(%rsp), %xmm4
	paddd	-104(%rsp), %xmm4
	movdqa	%xmm4, %xmm8
	psrld	$31, %xmm4
	pand	%xmm0, %xmm8
	paddd	%xmm8, %xmm4
	movups	%xmm4, 784(%rbx)
	movdqa	-104(%rsp), %xmm4
	pandn	%xmm0, %xmm4
	paddd	-120(%rsp), %xmm4
	movaps	%xmm7, -120(%rsp)
	movdqu	800(%rbx), %xmm7
	movdqa	%xmm4, %xmm8
	psrld	$31, %xmm4
	pand	%xmm0, %xmm8
	paddd	%xmm8, %xmm4
	movups	%xmm4, 656(%rbx)
	movdqa	%xmm3, %xmm4
	psrld	$31, %xmm3
	pand	%xmm0, %xmm4
	paddd	%xmm4, %xmm3
	movdqu	928(%rbx), %xmm4
	movups	%xmm3, 912(%rbx)
	movdqu	544(%rbx), %xmm3
	movaps	%xmm4, -104(%rsp)
	movdqa	%xmm7, %xmm4
	paddd	%xmm3, %xmm4
	pandn	%xmm0, %xmm3
	movdqa	%xmm4, %xmm8
	psrld	$31, %xmm4
	paddd	%xmm7, %xmm3
	pand	%xmm0, %xmm8
	paddd	%xmm8, %xmm4
	movups	%xmm4, 544(%rbx)
	movdqa	-104(%rsp), %xmm4
	paddd	-120(%rsp), %xmm4
	movdqa	%xmm4, %xmm8
	pand	%xmm0, %xmm4
	psrld	$31, %xmm8
	paddd	%xmm8, %xmm4
	movups	%xmm4, 800(%rbx)
	movdqa	-104(%rsp), %xmm4
	pandn	%xmm0, %xmm4
	paddd	-120(%rsp), %xmm4
	movdqa	%xmm4, %xmm8
	psrld	$31, %xmm4
	pand	%xmm0, %xmm8
	paddd	%xmm8, %xmm4
	movups	%xmm4, 672(%rbx)
	movdqa	%xmm3, %xmm4
	pand	%xmm0, %xmm3
	psrld	$31, %xmm4
	paddd	%xmm4, %xmm3
	movups	%xmm3, 928(%rbx)
	movdqu	560(%rbx), %xmm3
	movdqu	688(%rbx), %xmm7
	movdqu	944(%rbx), %xmm4
	movaps	%xmm7, -120(%rsp)
	movdqu	816(%rbx), %xmm7
	movaps	%xmm4, -104(%rsp)
	movdqa	%xmm3, %xmm4
	pandn	%xmm0, %xmm3
	paddd	%xmm7, %xmm4
	paddd	%xmm7, %xmm3
	movdqu	704(%rbx), %xmm7
	movdqa	%xmm4, %xmm8
	pand	%xmm0, %xmm4
	psrld	$31, %xmm8
	paddd	%xmm8, %xmm4
	movups	%xmm4, 560(%rbx)
	movdqa	-120(%rsp), %xmm4
	paddd	-104(%rsp), %xmm4
	movdqa	%xmm4, %xmm8
	pand	%xmm0, %xmm4
	psrld	$31, %xmm8
	paddd	%xmm8, %xmm4
	movups	%xmm4, 816(%rbx)
	movdqa	-104(%rsp), %xmm4
	pandn	%xmm0, %xmm4
	paddd	-120(%rsp), %xmm4
	movaps	%xmm7, -120(%rsp)
	movdqu	832(%rbx), %xmm7
	movdqa	%xmm4, %xmm8
	psrld	$31, %xmm4
	pand	%xmm0, %xmm8
	paddd	%xmm8, %xmm4
	movups	%xmm4, 688(%rbx)
	movdqa	%xmm3, %xmm4
	pand	%xmm0, %xmm3
	psrld	$31, %xmm4
	paddd	%xmm4, %xmm3
	movdqu	960(%rbx), %xmm4
	movups	%xmm3, 944(%rbx)
	movdqu	576(%rbx), %xmm3
	movaps	%xmm4, -104(%rsp)
	movdqa	%xmm3, %xmm4
	pandn	%xmm0, %xmm3
	paddd	%xmm7, %xmm4
	paddd	%xmm7, %xmm3
	movdqu	720(%rbx), %xmm7
	movdqa	%xmm4, %xmm8
	pand	%xmm0, %xmm4
	psrld	$31, %xmm8
	paddd	%xmm8, %xmm4
	movups	%xmm4, 576(%rbx)
	movdqa	-120(%rsp), %xmm4
	paddd	-104(%rsp), %xmm4
	movdqa	%xmm4, %xmm8
	psrld	$31, %xmm4
	pand	%xmm0, %xmm8
	paddd	%xmm8, %xmm4
	movups	%xmm4, 832(%rbx)
	movdqa	-104(%rsp), %xmm4
	pandn	%xmm0, %xmm4
	paddd	-120(%rsp), %xmm4
	movaps	%xmm7, -120(%rsp)
	movdqu	848(%rbx), %xmm7
	movdqa	%xmm4, %xmm8
	pand	%xmm0, %xmm4
	psrld	$31, %xmm8
	paddd	%xmm8, %xmm4
	movups	%xmm4, 704(%rbx)
	movdqa	%xmm3, %xmm4
	pand	%xmm0, %xmm3
	psrld	$31, %xmm4
	paddd	%xmm4, %xmm3
	movups	%xmm3, 960(%rbx)
	movdqu	592(%rbx), %xmm3
	movdqu	976(%rbx), %xmm4
	movaps	%xmm4, -104(%rsp)
	movdqa	%xmm3, %xmm4
	pandn	%xmm0, %xmm3
	paddd	%xmm7, %xmm4
	paddd	%xmm7, %xmm3
	movdqu	736(%rbx), %xmm7
	movdqa	%xmm4, %xmm8
	pand	%xmm0, %xmm4
	psrld	$31, %xmm8
	paddd	%xmm8, %xmm4
	movups	%xmm4, 592(%rbx)
	movdqa	-120(%rsp), %xmm4
	paddd	-104(%rsp), %xmm4
	movdqa	%xmm4, %xmm8
	pand	%xmm0, %xmm4
	psrld	$31, %xmm8
	paddd	%xmm8, %xmm4
	movups	%xmm4, 848(%rbx)
	movdqa	-104(%rsp), %xmm4
	pandn	%xmm0, %xmm4
	paddd	-120(%rsp), %xmm4
	movaps	%xmm7, -120(%rsp)
	movdqu	864(%rbx), %xmm7
	movdqa	%xmm4, %xmm8
	pand	%xmm0, %xmm4
	psrld	$31, %xmm8
	paddd	%xmm8, %xmm4
	movups	%xmm4, 720(%rbx)
	movdqa	%xmm3, %xmm4
	pand	%xmm0, %xmm3
	psrld	$31, %xmm4
	paddd	%xmm4, %xmm3
	movdqu	992(%rbx), %xmm4
	movups	%xmm3, 976(%rbx)
	movdqu	608(%rbx), %xmm3
	movaps	%xmm4, -104(%rsp)
	movdqa	%xmm7, %xmm4
	paddd	%xmm3, %xmm4
	pandn	%xmm0, %xmm3
	movdqa	%xmm4, %xmm8
	pand	%xmm0, %xmm4
	paddd	%xmm7, %xmm3
	movdqu	752(%rbx), %xmm7
	psrld	$31, %xmm8
	paddd	%xmm8, %xmm4
	movups	%xmm4, 608(%rbx)
	movdqa	-104(%rsp), %xmm4
	paddd	-120(%rsp), %xmm4
	movdqa	%xmm4, %xmm8
	psrld	$31, %xmm4
	pand	%xmm0, %xmm8
	paddd	%xmm8, %xmm4
	movups	%xmm4, 864(%rbx)
	movdqa	-104(%rsp), %xmm4
	pandn	%xmm0, %xmm4
	paddd	-120(%rsp), %xmm4
	movaps	%xmm7, -120(%rsp)
	movdqu	880(%rbx), %xmm7
	movdqa	%xmm4, %xmm8
	psrld	$31, %xmm4
	pand	%xmm0, %xmm8
	paddd	%xmm8, %xmm4
	movups	%xmm4, 736(%rbx)
	movdqa	%xmm3, %xmm4
	pand	%xmm0, %xmm3
	psrld	$31, %xmm4
	paddd	%xmm4, %xmm3
	movdqu	1008(%rbx), %xmm4
	movups	%xmm3, 992(%rbx)
	movdqu	624(%rbx), %xmm3
	movaps	%xmm4, -104(%rsp)
	movdqa	%xmm7, %xmm4
	paddd	%xmm3, %xmm4
	pandn	%xmm0, %xmm3
	movdqa	%xmm4, %xmm8
	psrld	$31, %xmm4
	paddd	%xmm7, %xmm3
	pand	%xmm0, %xmm8
	paddd	%xmm8, %xmm4
	movups	%xmm4, 624(%rbx)
	movdqa	-104(%rsp), %xmm4
	paddd	-120(%rsp), %xmm4
	movl	%esi, (%rbx)
	leal	(%rdi,%rdx), %esi
	movl	%esi, %r8d
	andl	$2147483647, %esi
	movdqa	%xmm4, %xmm8
	shrl	$31, %r8d
	psrld	$31, %xmm4
	pand	%xmm0, %xmm8
	addl	%r8d, %esi
	paddd	%xmm8, %xmm4
	leaq	omegas512(%rip), %r8
	movl	%esi, 512(%rbx)
	movl	%eax, %esi
	andl	$2147483647, %eax
	shrl	$31, %esi
	movups	%xmm4, 880(%rbx)
	movdqa	-104(%rsp), %xmm4
	addl	%eax, %esi
	movl	%edx, %eax
	notl	%eax
	pandn	%xmm0, %xmm4
	paddd	-120(%rsp), %xmm4
	andl	$2147483647, %eax
	addl	%edi, %eax
	movdqa	%xmm4, %xmm8
	movq	%r8, %rdi
	movl	%eax, %ecx
	andl	$2147483647, %eax
	pand	%xmm0, %xmm8
	shrl	$31, %ecx
	psrld	$31, %xmm4
	pand	%xmm3, %xmm0
	addl	%eax, %ecx
	psrld	$31, %xmm3
	paddd	%xmm8, %xmm4
	movl	%esi, %eax
	movl	%ecx, %edx
	sall	$15, %eax
	paddd	%xmm3, %xmm0
	movups	%xmm4, 752(%rbx)
	sall	$15, %edx
	shrl	$16, %esi
	andl	$2147450880, %eax
	movups	%xmm0, 1008(%rbx)
	shrl	$16, %ecx
	andl	$2147450880, %edx
	orl	%esi, %eax
	orl	%ecx, %edx
	leal	(%rdx,%rax), %ecx
	xorl	$2147483647, %eax
	addl	%edx, %eax
	movl	%ecx, %esi
	andl	$2147483647, %ecx
	movl	%eax, %edx
	andl	$2147483647, %eax
	shrl	$31, %esi
	shrl	$31, %edx
	addl	%esi, %ecx
	leaq	496(%r8), %rsi
	addl	%edx, %eax
	movl	%ecx, 256(%rbx)
	leaq	4(%rbx), %rcx
	movl	%eax, 768(%rbx)
	.p2align 4,,10
	.p2align 3
.L193:
	movl	256(%rcx), %r9d
	movl	768(%rcx), %r10d
	movl	(%rcx), %ebx
	movl	512(%rcx), %r11d
	leal	(%r9,%r10), %eax
	notl	%r10d
	movl	(%rdi), %ebp
	leal	(%rbx,%r11), %edx
	andl	$2147483647, %r10d
	notl	%r11d
	addl	%r9d, %r10d
	movl	%edx, %r9d
	andl	$2147483647, %edx
	andl	$2147483647, %r11d
	shrl	$31, %r9d
	addl	%ebx, %r11d
	leal	(%r9,%rdx), %r12d
	movq	%r12, %rbx
	imulq	%rbp, %rbx
	movq	%rbx, %rdx
	andl	$2147483647, %ebx
	shrq	$31, %rdx
	addl	%edx, %ebx
	movl	%eax, %edx
	andl	$2147483647, %eax
	shrl	$31, %edx
	addl	%edx, %eax
	movl	4(%rdi), %edx
	movl	%eax, %eax
	movq	%rax, %r9
	imulq	%rbp, %rax
	imulq	%rdx, %r9
	imulq	%r12, %rdx
	movq	%rax, %rbp
	andl	$2147483647, %eax
	movq	%r9, %r13
	andl	$2147483647, %r9d
	shrq	$31, %rbp
	shrq	$31, %r13
	movq	%rdx, %r12
	andl	$2147483647, %edx
	addl	%ebp, %eax
	addl	%r13d, %r9d
	shrq	$31, %r12
	addl	%r12d, %edx
	movl	%r9d, %ebp
	movl	%ebx, %r12d
	shrl	$31, %ebx
	andl	$2147483647, %r12d
	andl	$2147483647, %ebp
	shrl	$31, %r9d
	addl	%r12d, %ebp
	addl	%ebp, %ebx
	addl	%ebx, %r9d
	movl	%r9d, %ebx
	shrl	$31, %ebx
	andl	$2147483647, %r9d
	addl	%ebx, %r9d
	movl	%r9d, (%rcx)
	movl	%edx, %r9d
	andl	$2147483647, %edx
	shrl	$31, %r9d
	addl	%r9d, %edx
	movl	%eax, %r9d
	andl	$2147483647, %eax
	shrl	$31, %r9d
	notl	%edx
	andl	$2147483647, %edx
	addl	%r9d, %eax
	movq	%rsi, %r9
	addl	%edx, %eax
	movl	%eax, %edx
	andl	$2147483647, %eax
	shrl	$31, %edx
	addl	%edx, %eax
	movl	%eax, 512(%rcx)
	movl	%r11d, %eax
	andl	$2147483647, %r11d
	movl	4(%rsi), %ebx
	shrl	$31, %eax
	movl	(%rsi), %edx
	leal	(%rax,%r11), %ebp
	movq	%rbp, %r11
	imulq	%rbx, %r11
	movq	%r11, %rax
	andl	$2147483647, %r11d
	shrq	$31, %rax
	addl	%eax, %r11d
	movl	%r10d, %eax
	andl	$2147483647, %r10d
	shrl	$31, %eax
	addl	%r10d, %eax
	movl	%eax, %eax
	movq	%rax, %r10
	imulq	%rbx, %rax
	imulq	%rdx, %r10
	imulq	%rbp, %rdx
	movq	%rax, %rbx
	movq	%r10, %r12
	andl	$2147483647, %r10d
	movq	%rdx, %rbp
	shrq	$31, %r12
	andl	$2147483647, %edx
	shrq	$31, %rbp
	addl	%r12d, %r10d
	addl	%ebp, %edx
	shrq	$31, %rbx
	andl	$2147483647, %eax
	movl	%r11d, %ebp
	addl	%ebx, %eax
	movl	%r10d, %ebx
	andl	$2147483647, %ebp
	shrl	$31, %r11d
	shrl	$31, %r10d
	addq	$4, %rcx
	addq	$8, %rdi
	subq	$8, %rsi
	andl	$2147483647, %ebx
	addl	%ebp, %ebx
	addl	%ebx, %r11d
	addl	%r11d, %r10d
	movl	%r10d, %r11d
	andl	$2147483647, %r10d
	shrl	$31, %r11d
	addl	%r11d, %r10d
	movl	%r10d, 252(%rcx)
	movl	%edx, %r10d
	andl	$2147483647, %edx
	shrl	$31, %r10d
	addl	%r10d, %edx
	movl	%eax, %r10d
	andl	$2147483647, %eax
	notl	%edx
	shrl	$31, %r10d
	andl	$2147483647, %edx
	addl	%r10d, %eax
	addl	%edx, %eax
	movl	%eax, %edx
	andl	$2147483647, %eax
	shrl	$31, %edx
	addl	%edx, %eax
	movl	%eax, 764(%rcx)
	cmpq	%r8, %r9
	jne	.L193
	addq	$64, %rsp
	.cfi_remember_state
	.cfi_def_cfa_offset 56
	popq	%rbx
	.cfi_def_cfa_offset 48
	popq	%rbp
	.cfi_def_cfa_offset 40
	popq	%r12
	.cfi_def_cfa_offset 32
	popq	%r13
	.cfi_def_cfa_offset 24
	popq	%r14
	.cfi_def_cfa_offset 16
	popq	%r15
	.cfi_def_cfa_offset 8
	ret
	.p2align 4,,10
	.p2align 3
.L198:
	.cfi_restore_state
	movq	%r14, %r9
	.p2align 4,,10
	.p2align 3
.L188:
	movl	(%r9), %r11d
	movl	(%r9,%r15,4), %r8d
	movl	(%r9,%r13,4), %r10d
	movl	(%r9,%rsi,4), %edi
	leal	(%r8,%r11), %ebx
	notl	%r8d
	movl	%ebx, %ebp
	shrl	$31, %ebx
	andl	$2147483647, %r8d
	andl	$2147483647, %ebp
	addl	%r11d, %r8d
	addl	%ebx, %ebp
	leal	(%rdi,%r10), %ebx
	notl	%edi
	andl	$2147483647, %edi
	movl	%ebp, (%r9)
	movl	%ebx, %ebp
	shrl	$31, %ebx
	addl	%r10d, %edi
	movl	%r8d, %r10d
	shrl	$31, %r8d
	andl	$2147483647, %ebp
	andl	$2147483647, %r10d
	addl	%ebx, %ebp
	addl	%r10d, %r8d
	movl	%ebp, (%r9,%r15,4)
	movl	%r8d, %r8d
	movq	%r8, %r10
	imulq	%rax, %r10
	movl	%r10d, %r11d
	shrq	$31, %r10
	andl	$2147483647, %r11d
	addl	%r10d, %r11d
	movl	%edi, %r10d
	shrl	$31, %edi
	andl	$2147483647, %r10d
	movl	%r11d, %ebp
	shrl	$31, %r11d
	addl	%r10d, %edi
	andl	$2147483647, %ebp
	movl	%edi, %edi
	movq	%rdi, %rbx
	imulq	%rdx, %rbx
	movl	%ebx, %r10d
	shrq	$31, %rbx
	andl	$2147483647, %r10d
	addl	%ebx, %r10d
	movl	%r10d, %ebx
	shrl	$31, %r10d
	andl	$2147483647, %ebx
	addl	%ebp, %ebx
	addl	%ebx, %r11d
	addl	%r11d, %r10d
	movl	%r10d, %r11d
	andl	$2147483647, %r11d
	shrl	$31, %r10d
	imulq	%rdx, %r8
	addl	%r11d, %r10d
	imulq	%rax, %rdi
	movl	%r10d, (%r9,%r13,4)
	movl	%r8d, %r10d
	shrq	$31, %r8
	andl	$2147483647, %r10d
	addl	%r8d, %r10d
	movl	%edi, %r8d
	shrq	$31, %rdi
	andl	$2147483647, %r8d
	movl	%r10d, %r11d
	addl	%edi, %r8d
	movl	%r10d, %edi
	andl	$2147483647, %r11d
	shrl	$31, %edi
	movl	%r8d, %r10d
	shrl	$31, %r8d
	addl	%r11d, %edi
	andl	$2147483647, %r10d
	notl	%edi
	addl	%r10d, %r8d
	andl	$2147483647, %edi
	addl	%r8d, %edi
	movl	%edi, %r8d
	shrl	$31, %edi
	andl	$2147483647, %r8d
	addl	%r8d, %edi
	movl	%edi, (%r9,%rsi,4)
	addq	$4, %r9
	cmpq	%r9, %rcx
	jne	.L188
	jmp	.L192
	.p2align 4,,10
	.p2align 3
.L144:
	movslq	-68(%rsp), %rsi
	movq	%r14, %rcx
	leaq	(%r14,%r13,4), %r11
	.p2align 4,,10
	.p2align 3
.L152:
	movl	(%rcx), %r9d
	movl	(%rcx,%r15,4), %edx
	movl	(%rcx,%r13,4), %edi
	movl	(%rcx,%rsi,4), %eax
	leal	(%r9,%rdx), %r8d
	notl	%edx
	movl	%r8d, %r10d
	andl	$2147483647, %r8d
	andl	$2147483647, %edx
	shrl	$31, %r10d
	addl	%r9d, %edx
	addl	%r10d, %r8d
	movl	%r8d, (%rcx)
	leal	(%rdi,%rax), %r8d
	notl	%eax
	movl	%r8d, %r10d
	andl	$2147483647, %r8d
	andl	$2147483647, %eax
	shrl	$31, %r10d
	addl	%edi, %eax
	addl	%r10d, %r8d
	movl	%r8d, (%rcx,%r15,4)
	movl	%edx, %r8d
	andl	$2147483647, %edx
	shrl	$31, %r8d
	addl	%r8d, %edx
	movl	%edx, (%rcx,%r13,4)
	movl	%eax, %edx
	andl	$2147483647, %eax
	shrl	$31, %edx
	addl	%edx, %eax
	movl	%eax, (%rcx,%rsi,4)
	addq	$4, %rcx
	cmpq	%rcx, %r11
	jne	.L152
	cmpl	$1, -72(%rsp)
	movslq	-36(%rsp), %r12
	jne	.L149
.L154:
	leaq	0(%r13,%r12), %rax
	leaq	(%r14,%r12,4), %rcx
	leaq	(%r14,%rax,4), %r11
	.p2align 4,,10
	.p2align 3
.L162:
	movl	(%rcx), %eax
	movl	(%rcx,%r15,4), %edi
	movl	(%rcx,%r13,4), %r9d
	movl	(%rcx,%rsi,4), %edx
	leal	(%rax,%rdi), %r8d
	notl	%eax
	movl	%r8d, %r10d
	andl	$2147483647, %r8d
	andl	$2147483647, %eax
	shrl	$31, %r10d
	addl	%edi, %eax
	addl	%r10d, %r8d
	movl	%r8d, (%rcx)
	leal	(%r9,%rdx), %r8d
	notl	%edx
	movl	%r8d, %r10d
	andl	$2147483647, %edx
	andl	$2147483647, %r8d
	shrl	$31, %r10d
	addl	%r9d, %edx
	addl	%r10d, %r8d
	movl	%r8d, (%rcx,%r15,4)
	movl	%edx, %r8d
	andl	$2147483647, %edx
	shrl	$31, %r8d
	addl	%r8d, %edx
	movl	%edx, (%rcx,%r13,4)
	movl	%eax, %edx
	andl	$2147483647, %eax
	shrl	$31, %edx
	addl	%edx, %eax
	movl	%eax, (%rcx,%rsi,4)
	addq	$4, %rcx
	cmpq	%r11, %rcx
	jne	.L162
	cmpl	$1, -72(%rsp)
	jne	.L159
.L164:
	movslq	-120(%rsp), %rax
	leaq	(%r14,%rax,4), %rcx
	addq	%r13, %rax
	leaq	(%r14,%rax,4), %r11
	.p2align 4,,10
	.p2align 3
.L172:
	movl	(%rcx), %r9d
	movl	(%rcx,%r15,4), %edx
	movl	(%rcx,%r13,4), %edi
	movl	(%rcx,%rsi,4), %eax
	leal	(%r9,%rdx), %r8d
	notl	%edx
	movl	%r8d, %r10d
	andl	$2147483647, %r8d
	andl	$2147483647, %edx
	shrl	$31, %r10d
	addl	%r9d, %edx
	addl	%r10d, %r8d
	movl	%r8d, (%rcx)
	leal	(%rdi,%rax), %r8d
	notl	%eax
	movl	%r8d, %r10d
	andl	$2147483647, %r8d
	andl	$2147483647, %eax
	shrl	$31, %r10d
	addl	%edi, %eax
	addl	%r10d, %r8d
	movl	%eax, %edi
	andl	$2147483647, %eax
	movl	%r8d, (%rcx,%r15,4)
	movl	%edx, %r8d
	andl	$2147483647, %edx
	shrl	$31, %edi
	shrl	$31, %r8d
	addl	%eax, %edi
	addl	%r8d, %edx
	movl	%edx, %eax
	shrl	$16, %edx
	sall	$15, %eax
	andl	$2147450880, %eax
	orl	%edx, %eax
	movl	%edi, %edx
	shrl	$16, %edi
	sall	$15, %edx
	andl	$2147450880, %edx
	orl	%edi, %edx
	leal	(%rax,%rdx), %edi
	xorl	$2147483647, %eax
	addl	%edx, %eax
	movl	%edi, %r8d
	andl	$2147483647, %edi
	movl	%eax, %edx
	shrl	$31, %r8d
	andl	$2147483647, %eax
	shrl	$31, %edx
	addl	%r8d, %edi
	addl	%edx, %eax
	movl	%edi, (%rcx,%r13,4)
	movl	%eax, (%rcx,%rsi,4)
	addq	$4, %rcx
	cmpq	%r11, %rcx
	jne	.L172
	movl	-120(%rsp), %eax
	movl	-36(%rsp), %ecx
	addl	%ecx, %eax
	cmpl	$1, -72(%rsp)
	movl	%eax, -64(%rsp)
	jne	.L169
.L174:
	movslq	-64(%rsp), %rax
	leaq	(%r14,%rax,4), %rcx
	addq	%r13, %rax
	leaq	(%r14,%rax,4), %r11
	.p2align 4,,10
	.p2align 3
.L180:
	movl	(%rcx), %edx
	movl	(%rcx,%r15,4), %r9d
	movl	(%rcx,%r13,4), %edi
	movl	(%rcx,%rsi,4), %eax
	leal	(%r9,%rdx), %r8d
	notl	%edx
	movl	%r8d, %r10d
	shrl	$31, %r8d
	andl	$2147483647, %edx
	andl	$2147483647, %r10d
	addl	%r9d, %edx
	addl	%r10d, %r8d
	movl	%r8d, (%rcx)
	leal	(%rax,%rdi), %r8d
	notl	%eax
	movl	%r8d, %r10d
	andl	$2147483647, %eax
	shrl	$31, %r8d
	andl	$2147483647, %r10d
	addl	%edi, %eax
	addl	%r10d, %r8d
	movl	%eax, %edi
	shrl	$31, %eax
	movl	%r8d, (%rcx,%r15,4)
	movl	%edx, %r8d
	andl	$2147483647, %edi
	shrl	$31, %edx
	andl	$2147483647, %r8d
	addl	%eax, %edi
	addl	%edx, %r8d
	movl	%edi, %eax
	shrl	$16, %edi
	movl	%r8d, %edx
	sall	$15, %eax
	sall	$15, %edx
	shrl	$16, %r8d
	andl	$2147450880, %eax
	andl	$2147450880, %edx
	orl	%edi, %eax
	orl	%r8d, %edx
	leal	(%rax,%rdx), %edi
	xorl	$2147483647, %eax
	addl	%edx, %eax
	movl	%edi, %r8d
	shrl	$31, %edi
	movl	%eax, %edx
	andl	$2147483647, %r8d
	shrl	$31, %eax
	andl	$2147483647, %edx
	addl	%r8d, %edi
	addl	%edx, %eax
	movl	%edi, (%rcx,%r13,4)
	movl	%eax, (%rcx,%rsi,4)
	addq	$4, %rcx
	cmpq	%r11, %rcx
	jne	.L180
	jmp	.L181
.L302:
	movslq	-68(%rsp), %rsi
	jmp	.L154
.L304:
	movslq	-68(%rsp), %rsi
	jmp	.L174
.L303:
	movslq	-68(%rsp), %rsi
	jmp	.L164
.L166:
	movl	-72(%rsp), %esi
	movl	%esi, %eax
	andl	$-4, %eax
	movl	%eax, %edx
	cmpl	%eax, %esi
	je	.L171
	movl	%esi, %r8d
	subl	%eax, %r8d
	cmpl	$1, %r8d
	je	.L170
.L165:
	leaq	0(%rbp,%rdx), %rcx
	movq	-80(%rsp), %rbx
	movq	.LC4(%rip), %xmm7
	leaq	(%r14,%rcx,4), %rdi
	movq	-56(%rsp), %rcx
	movq	(%rdi), %xmm6
	addq	%rdx, %rbx
	addq	%rdx, %rcx
	leaq	(%r14,%rcx,4), %rsi
	movq	-64(%rsp), %rcx
	movdqa	%xmm6, %xmm3
	movq	(%rsi), %xmm4
	addq	%rdx, %rcx
	leaq	(%r14,%rbx,4), %rdx
	leaq	(%r14,%rcx,4), %rcx
	movq	(%rdx), %xmm1
	movq	(%rcx), %xmm2
	paddd	%xmm2, %xmm3
	pandn	%xmm14, %xmm2
	paddd	%xmm6, %xmm2
	movdqa	%xmm3, %xmm5
	pand	%xmm14, %xmm3
	psrld	$31, %xmm5
	paddd	%xmm5, %xmm3
	movq	%xmm3, (%rdi)
	movdqa	%xmm4, %xmm3
	paddd	%xmm1, %xmm3
	pandn	%xmm14, %xmm1
	paddd	%xmm4, %xmm1
	movdqa	%xmm3, %xmm5
	pand	%xmm14, %xmm3
	psrld	$31, %xmm5
	paddd	%xmm5, %xmm3
	movq	%xmm3, (%rcx)
	movdqa	%xmm2, %xmm3
	pand	%xmm14, %xmm2
	psrld	$31, %xmm3
	paddd	%xmm3, %xmm2
	movdqa	%xmm1, %xmm3
	pand	%xmm14, %xmm1
	psrld	$31, %xmm3
	paddd	%xmm1, %xmm3
	movdqa	%xmm2, %xmm1
	psrld	$16, %xmm2
	pslld	$15, %xmm1
	pand	%xmm7, %xmm1
	por	%xmm2, %xmm1
	movdqa	%xmm3, %xmm2
	pslld	$15, %xmm2
	psrld	$16, %xmm3
	pand	%xmm7, %xmm2
	por	%xmm3, %xmm2
	movdqa	%xmm1, %xmm3
	pxor	%xmm14, %xmm1
	paddd	%xmm2, %xmm3
	paddd	%xmm2, %xmm1
	movdqa	%xmm3, %xmm4
	movdqa	%xmm1, %xmm2
	pand	%xmm14, %xmm3
	psrld	$31, %xmm4
	psrld	$31, %xmm2
	pand	%xmm14, %xmm1
	paddd	%xmm4, %xmm3
	paddd	%xmm2, %xmm1
	movq	%xmm3, (%rsi)
	movq	%xmm1, (%rdx)
	testb	$1, %r8b
	je	.L171
	andl	$-2, %r8d
	addl	%r8d, %eax
.L170:
	movl	-120(%rsp), %esi
	addl	%esi, %eax
	movl	-72(%rsp), %esi
	movslq	%eax, %rdx
	leaq	(%r14,%rdx,4), %rbx
	leal	(%rsi,%rax), %edx
	movl	-40(%rsp), %esi
	movslq	%edx, %rdx
	movl	(%rbx), %r10d
	leaq	(%r14,%rdx,4), %r8
	leal	(%rsi,%rax), %edx
	movl	-68(%rsp), %esi
	movslq	%edx, %rdx
	movl	(%r8), %r9d
	leaq	(%r14,%rdx,4), %r11
	addl	%esi, %eax
	movl	(%r11), %edx
	cltq
	leaq	(%r14,%rax,4), %rdi
	leal	(%r10,%rdx), %ecx
	movl	(%rdi), %eax
	notl	%edx
	movl	%ecx, %esi
	andl	$2147483647, %ecx
	andl	$2147483647, %edx
	shrl	$31, %esi
	addl	%r10d, %edx
	addl	%ecx, %esi
	leal	(%r9,%rax), %ecx
	notl	%eax
	movl	%esi, (%rbx)
	movl	%ecx, %esi
	andl	$2147483647, %ecx
	andl	$2147483647, %eax
	shrl	$31, %esi
	addl	%ecx, %esi
	movl	%edx, %ecx
	andl	$2147483647, %edx
	shrl	$31, %ecx
	movl	%esi, (%r11)
	leal	(%rax,%r9), %esi
	addl	%edx, %ecx
	movl	%esi, %eax
	andl	$2147483647, %esi
	movl	%ecx, %edx
	shrl	$31, %eax
	sall	$15, %edx
	addl	%eax, %esi
	shrl	$16, %ecx
	movl	%edx, %eax
	movl	%ecx, %edx
	andl	$2147450880, %eax
	orl	%eax, %edx
	movl	%esi, %eax
	shrl	$16, %esi
	sall	$15, %eax
	andl	$2147450880, %eax
	orl	%esi, %eax
	leal	(%rdx,%rax), %ecx
	xorl	$2147483647, %edx
	addl	%eax, %edx
	movl	%ecx, %esi
	andl	$2147483647, %ecx
	movl	%edx, %eax
	shrl	$31, %esi
	andl	$2147483647, %edx
	shrl	$31, %eax
	addl	%ecx, %esi
	addl	%edx, %eax
	movl	%esi, (%r8)
	movl	%eax, (%rdi)
	jmp	.L171
.L156:
	movl	-72(%rsp), %esi
	movl	%esi, %edx
	andl	$-4, %edx
	movl	%edx, %eax
	cmpl	%edx, %esi
	je	.L159
	subl	%edx, %esi
	movl	%esi, %ecx
	cmpl	$1, %esi
	je	.L160
.L155:
	leaq	(%r12,%rax), %rsi
	leaq	0(%rbp,%rax), %r8
	movq	.LC1(%rip), %xmm3
	leaq	(%r14,%rsi,4), %rdi
	leaq	(%r14,%r8,4), %r8
	movq	-64(%rsp), %rsi
	movq	(%rdi), %xmm1
	movq	(%r8), %xmm5
	addq	%rax, %rsi
	addq	%rbx, %rax
	movdqa	%xmm1, %xmm4
	leaq	(%r14,%rsi,4), %rsi
	leaq	(%r14,%rax,4), %rax
	paddd	%xmm5, %xmm4
	movq	(%rsi), %xmm7
	movq	(%rax), %xmm2
	pandn	%xmm3, %xmm1
	paddd	%xmm5, %xmm1
	movdqa	%xmm4, %xmm6
	pand	%xmm3, %xmm4
	psrld	$31, %xmm6
	paddd	%xmm6, %xmm4
	movq	%xmm4, (%rdi)
	movdqa	%xmm7, %xmm4
	paddd	%xmm2, %xmm4
	pandn	%xmm3, %xmm2
	paddd	%xmm7, %xmm2
	movdqa	%xmm4, %xmm6
	pand	%xmm3, %xmm4
	psrld	$31, %xmm6
	paddd	%xmm6, %xmm4
	movq	%xmm4, (%r8)
	movdqa	%xmm2, %xmm4
	pand	%xmm3, %xmm2
	psrld	$31, %xmm4
	paddd	%xmm4, %xmm2
	movq	%xmm2, (%rsi)
	movdqa	%xmm1, %xmm2
	pand	%xmm3, %xmm1
	psrld	$31, %xmm2
	paddd	%xmm2, %xmm1
	movq	%xmm1, (%rax)
	testb	$1, %cl
	je	.L159
	andl	$-2, %ecx
	addl	%ecx, %edx
.L160:
	movl	-36(%rsp), %eax
	movl	-72(%rsp), %esi
	addl	%eax, %edx
	leal	(%rsi,%rdx), %ecx
	movl	-40(%rsp), %esi
	movslq	%edx, %rax
	movslq	%ecx, %rcx
	leaq	(%r14,%rax,4), %rbx
	leaq	(%r14,%rcx,4), %r9
	leal	(%rsi,%rdx), %ecx
	movl	-68(%rsp), %esi
	movl	(%rbx), %eax
	movslq	%ecx, %rcx
	movl	(%r9), %r10d
	leaq	(%r14,%rcx,4), %r11
	addl	%esi, %edx
	movl	(%r11), %r8d
	movslq	%edx, %rdx
	leaq	(%r14,%rdx,4), %rdi
	leal	(%rax,%r8), %ecx
	movl	(%rdi), %edx
	notl	%eax
	movl	%ecx, %esi
	andl	$2147483647, %ecx
	andl	$2147483647, %eax
	shrl	$31, %esi
	addl	%r8d, %eax
	addl	%ecx, %esi
	leal	(%r10,%rdx), %ecx
	notl	%edx
	movl	%esi, (%rbx)
	andl	$2147483647, %edx
	movl	%ecx, %esi
	andl	$2147483647, %ecx
	addl	%r10d, %edx
	shrl	$31, %esi
	addl	%ecx, %esi
	movl	%edx, %ecx
	andl	$2147483647, %edx
	shrl	$31, %ecx
	movl	%esi, (%r11)
	addl	%edx, %ecx
	movl	%eax, %edx
	andl	$2147483647, %eax
	shrl	$31, %edx
	movl	%ecx, (%r9)
	addl	%edx, %eax
	movl	%eax, (%rdi)
	jmp	.L159
.L176:
	movl	-72(%rsp), %esi
	movl	%esi, %eax
	andl	$-4, %eax
	movl	%eax, %edx
	cmpl	%eax, %esi
	je	.L181
	movl	%esi, %r8d
	subl	%eax, %r8d
	cmpl	$1, %r8d
	je	.L178
.L175:
	movq	.LC4(%rip), %xmm7
	leaq	(%rbx,%rdx), %rcx
	leaq	(%r14,%rcx,4), %rdi
	movq	-80(%rsp), %rcx
	movq	(%rdi), %xmm2
	addq	%rdx, %rcx
	leaq	(%r14,%rcx,4), %rsi
	movq	-56(%rsp), %rcx
	movdqa	%xmm2, %xmm3
	pandn	%xmm14, %xmm2
	movq	(%rsi), %xmm5
	addq	%rdx, %rcx
	addq	%rbp, %rdx
	leaq	(%r14,%rcx,4), %rcx
	leaq	(%r14,%rdx,4), %rdx
	movq	(%rcx), %xmm6
	movq	(%rdx), %xmm1
	paddd	%xmm6, %xmm3
	paddd	%xmm6, %xmm2
	movdqa	%xmm3, %xmm4
	pand	%xmm14, %xmm3
	psrld	$31, %xmm4
	paddd	%xmm4, %xmm3
	movq	%xmm3, (%rdi)
	movdqa	%xmm5, %xmm3
	paddd	%xmm1, %xmm3
	pandn	%xmm14, %xmm1
	paddd	%xmm5, %xmm1
	movdqa	%xmm3, %xmm4
	pand	%xmm14, %xmm3
	psrld	$31, %xmm4
	paddd	%xmm4, %xmm3
	movdqa	%xmm2, %xmm4
	pand	%xmm14, %xmm2
	psrld	$31, %xmm4
	paddd	%xmm2, %xmm4
	movq	%xmm3, (%rcx)
	movdqa	%xmm1, %xmm3
	pand	%xmm14, %xmm1
	psrld	$31, %xmm3
	paddd	%xmm1, %xmm3
	movdqa	%xmm4, %xmm2
	pslld	$15, %xmm2
	psrld	$16, %xmm4
	pand	%xmm7, %xmm2
	movdqa	%xmm3, %xmm1
	psrld	$16, %xmm3
	por	%xmm4, %xmm2
	pslld	$15, %xmm1
	pand	%xmm7, %xmm1
	por	%xmm3, %xmm1
	movdqa	%xmm2, %xmm3
	paddd	%xmm1, %xmm3
	pxor	%xmm14, %xmm1
	paddd	%xmm2, %xmm1
	movdqa	%xmm3, %xmm4
	pand	%xmm14, %xmm3
	movdqa	%xmm1, %xmm2
	psrld	$31, %xmm4
	pand	%xmm14, %xmm1
	psrld	$31, %xmm2
	paddd	%xmm4, %xmm3
	paddd	%xmm2, %xmm1
	movq	%xmm3, (%rsi)
	movq	%xmm1, (%rdx)
	testb	$1, %r8b
	je	.L181
	andl	$-2, %r8d
	addl	%r8d, %eax
.L178:
	movl	-64(%rsp), %esi
	addl	%esi, %eax
	movl	-72(%rsp), %esi
	movslq	%eax, %rdx
	leal	(%rsi,%rax), %ecx
	movl	-40(%rsp), %esi
	leaq	(%r14,%rdx,4), %rbx
	movslq	%ecx, %rcx
	movl	(%rbx), %edx
	leaq	(%r14,%rcx,4), %r8
	leal	(%rsi,%rax), %ecx
	movl	-68(%rsp), %esi
	movslq	%ecx, %rcx
	movl	(%r8), %r9d
	leaq	(%r14,%rcx,4), %r11
	addl	%esi, %eax
	movl	(%r11), %r10d
	cltq
	leaq	(%r14,%rax,4), %rdi
	leal	(%rdx,%r10), %ecx
	movl	(%rdi), %eax
	notl	%edx
	movl	%ecx, %esi
	andl	$2147483647, %ecx
	andl	$2147483647, %edx
	shrl	$31, %esi
	addl	%r10d, %edx
	addl	%ecx, %esi
	leal	(%r9,%rax), %ecx
	notl	%eax
	movl	%esi, (%rbx)
	movl	%ecx, %esi
	andl	$2147483647, %ecx
	andl	$2147483647, %eax
	shrl	$31, %esi
	addl	%ecx, %esi
	movl	%edx, %ecx
	andl	$2147483647, %edx
	shrl	$31, %ecx
	movl	%esi, (%r11)
	leal	(%rax,%r9), %esi
	addl	%edx, %ecx
	movl	%esi, %eax
	andl	$2147483647, %esi
	movl	%ecx, %edx
	shrl	$31, %eax
	sall	$15, %edx
	addl	%eax, %esi
	shrl	$16, %ecx
	movl	%edx, %eax
	movl	%ecx, %edx
	andl	$2147450880, %eax
	orl	%eax, %edx
	movl	%esi, %eax
	shrl	$16, %esi
	sall	$15, %eax
	andl	$2147450880, %eax
	orl	%esi, %eax
	leal	(%rdx,%rax), %ecx
	xorl	$2147483647, %eax
	addl	%edx, %eax
	movl	%ecx, %esi
	andl	$2147483647, %ecx
	movl	%eax, %edx
	shrl	$31, %esi
	andl	$2147483647, %eax
	shrl	$31, %edx
	addl	%ecx, %esi
	addl	%edx, %eax
	movl	%esi, (%r8)
	movl	%eax, (%rdi)
	jmp	.L181
.L146:
	movl	-72(%rsp), %esi
	movl	%esi, %eax
	andl	$-4, %eax
	movl	%eax, %edx
	cmpl	%eax, %esi
	je	.L149
	subl	%eax, %esi
	movl	%esi, %ecx
	cmpl	$1, %esi
	je	.L150
.L145:
	leaq	(%r14,%rdx,4), %r8
	leaq	(%r15,%rdx), %rdi
	addq	%rdx, %r12
	movq	.LC1(%rip), %xmm3
	movq	(%r8), %xmm7
	leaq	(%r14,%rdi,4), %rdi
	leaq	0(%r13,%rdx), %rsi
	movq	(%rdi), %xmm2
	leaq	(%r14,%rsi,4), %rsi
	leaq	(%r14,%r12,4), %rdx
	movdqa	%xmm7, %xmm4
	movq	(%rsi), %xmm5
	movq	(%rdx), %xmm1
	paddd	%xmm2, %xmm4
	pandn	%xmm3, %xmm2
	paddd	%xmm7, %xmm2
	movdqa	%xmm4, %xmm6
	pand	%xmm3, %xmm4
	psrld	$31, %xmm6
	paddd	%xmm6, %xmm4
	movq	%xmm4, (%r8)
	movdqa	%xmm5, %xmm4
	paddd	%xmm1, %xmm4
	pandn	%xmm3, %xmm1
	paddd	%xmm5, %xmm1
	movdqa	%xmm4, %xmm6
	pand	%xmm3, %xmm4
	psrld	$31, %xmm6
	paddd	%xmm6, %xmm4
	movq	%xmm4, (%rdi)
	movdqa	%xmm2, %xmm4
	pand	%xmm3, %xmm2
	psrld	$31, %xmm4
	paddd	%xmm4, %xmm2
	movq	%xmm2, (%rsi)
	movdqa	%xmm1, %xmm2
	pand	%xmm3, %xmm1
	psrld	$31, %xmm2
	paddd	%xmm2, %xmm1
	movq	%xmm1, (%rdx)
	testb	$1, %cl
	je	.L149
	andl	$-2, %ecx
	addl	%ecx, %eax
.L150:
	movl	-72(%rsp), %esi
	movslq	%eax, %rdx
	leaq	(%r14,%rdx,4), %rbx
	leal	(%rsi,%rax), %edx
	movl	-40(%rsp), %esi
	movl	(%rbx), %r10d
	movslq	%edx, %rdx
	leaq	(%r14,%rdx,4), %r9
	leal	(%rsi,%rax), %edx
	movl	-68(%rsp), %esi
	movslq	%edx, %rdx
	movl	(%r9), %r8d
	leaq	(%r14,%rdx,4), %r11
	addl	%esi, %eax
	movl	(%r11), %edx
	cltq
	leaq	(%r14,%rax,4), %rdi
	leal	(%r10,%rdx), %ecx
	movl	(%rdi), %eax
	notl	%edx
	movl	%ecx, %esi
	andl	$2147483647, %ecx
	andl	$2147483647, %edx
	shrl	$31, %esi
	addl	%r10d, %edx
	addl	%ecx, %esi
	leal	(%r8,%rax), %ecx
	notl	%eax
	movl	%esi, (%rbx)
	movl	%ecx, %esi
	andl	$2147483647, %ecx
	andl	$2147483647, %eax
	shrl	$31, %esi
	addl	%r8d, %eax
	addl	%ecx, %esi
	movl	%edx, %ecx
	andl	$2147483647, %edx
	shrl	$31, %ecx
	movl	%esi, (%r11)
	addl	%edx, %ecx
	movl	%eax, %edx
	andl	$2147483647, %eax
	shrl	$31, %edx
	movl	%ecx, (%r9)
	addl	%edx, %eax
	movl	%eax, (%rdi)
	jmp	.L149
.L194:
	xorl	%edx, %edx
	xorl	%eax, %eax
	jmp	.L145
.L195:
	xorl	%eax, %eax
	xorl	%edx, %edx
	jmp	.L155
.L196:
	xorl	%edx, %edx
	xorl	%eax, %eax
	jmp	.L165
.L197:
	xorl	%edx, %edx
	xorl	%eax, %eax
	jmp	.L175
	.cfi_endproc
.LFE50:
	.size	mrsn_invntt_256, .-mrsn_invntt_256
	.p2align 4
	.globl	mrsn_mulc_256
	.type	mrsn_mulc_256, @function
mrsn_mulc_256:
.LFB51:
	.cfi_startproc
	endbr64
	movq	%rdi, %rax
	movq	%rdx, %r8
	subq	%rdx, %rax
	subq	$4, %rax
	cmpq	$24, %rax
	jbe	.L309
	movq	%rdi, %rax
	subq	%rsi, %rax
	subq	$4, %rax
	cmpq	$24, %rax
	jbe	.L309
	movdqa	.LC0(%rip), %xmm0
	xorl	%eax, %eax
	.p2align 4,,10
	.p2align 3
.L307:
	movdqu	(%rsi,%rax), %xmm7
	movdqu	16(%rsi,%rax), %xmm5
	movdqu	(%r8,%rax), %xmm3
	movdqa	%xmm7, %xmm9
	shufps	$221, %xmm5, %xmm7
	shufps	$136, %xmm5, %xmm9
	movdqu	16(%r8,%rax), %xmm5
	movdqa	%xmm3, %xmm6
	movdqa	%xmm9, %xmm10
	punpckldq	%xmm9, %xmm10
	punpckhdq	%xmm9, %xmm9
	shufps	$136, %xmm5, %xmm6
	movdqa	%xmm6, %xmm1
	shufps	$221, %xmm5, %xmm3
	punpckldq	%xmm6, %xmm1
	punpckhdq	%xmm6, %xmm6
	movdqa	%xmm1, %xmm2
	movdqa	%xmm6, %xmm4
	pmuludq	%xmm9, %xmm4
	pmuludq	%xmm10, %xmm2
	movdqa	%xmm4, %xmm8
	movdqa	%xmm2, %xmm5
	shufps	$136, %xmm4, %xmm2
	pand	%xmm0, %xmm2
	psrlq	$31, %xmm8
	psrlq	$31, %xmm5
	shufps	$136, %xmm8, %xmm5
	paddd	%xmm2, %xmm5
	movdqa	%xmm3, %xmm2
	movdqa	%xmm7, %xmm8
	punpckldq	%xmm3, %xmm2
	punpckhdq	%xmm3, %xmm3
	punpckldq	%xmm7, %xmm8
	movdqa	%xmm2, %xmm4
	movdqa	%xmm3, %xmm12
	punpckhdq	%xmm7, %xmm7
	pmuludq	%xmm9, %xmm3
	pmuludq	%xmm10, %xmm2
	pmuludq	%xmm7, %xmm12
	pmuludq	%xmm8, %xmm4
	pmuludq	%xmm8, %xmm1
	movdqa	%xmm2, %xmm9
	movdqa	%xmm3, %xmm10
	shufps	$136, %xmm3, %xmm2
	pand	%xmm0, %xmm2
	psrlq	$31, %xmm9
	psrlq	$31, %xmm10
	movdqa	%xmm4, %xmm11
	shufps	$136, %xmm12, %xmm4
	shufps	$136, %xmm10, %xmm9
	movdqa	%xmm9, %xmm3
	movdqa	%xmm12, %xmm13
	pand	%xmm0, %xmm4
	psrlq	$31, %xmm11
	psrlq	$31, %xmm13
	paddd	%xmm2, %xmm3
	movdqa	%xmm6, %xmm2
	movdqa	%xmm1, %xmm6
	shufps	$136, %xmm13, %xmm11
	paddd	%xmm11, %xmm4
	pmuludq	%xmm7, %xmm2
	psrlq	$31, %xmm6
	movdqa	%xmm2, %xmm7
	shufps	$136, %xmm2, %xmm1
	pand	%xmm0, %xmm1
	movdqa	%xmm5, %xmm2
	psrlq	$31, %xmm7
	psrld	$31, %xmm2
	pand	%xmm0, %xmm5
	shufps	$136, %xmm7, %xmm6
	paddd	%xmm1, %xmm6
	movdqa	%xmm4, %xmm1
	pand	%xmm0, %xmm4
	psrld	$31, %xmm1
	paddd	%xmm5, %xmm2
	paddd	%xmm4, %xmm1
	pandn	%xmm0, %xmm1
	paddd	%xmm2, %xmm1
	movdqa	%xmm1, %xmm2
	pand	%xmm0, %xmm1
	psrld	$31, %xmm2
	paddd	%xmm1, %xmm2
	movdqa	%xmm3, %xmm1
	pand	%xmm0, %xmm3
	psrld	$31, %xmm1
	paddd	%xmm3, %xmm1
	movdqa	%xmm6, %xmm3
	pand	%xmm0, %xmm6
	psrld	$31, %xmm3
	paddd	%xmm6, %xmm3
	paddd	%xmm3, %xmm1
	movdqa	%xmm1, %xmm3
	pand	%xmm0, %xmm1
	psrld	$31, %xmm3
	paddd	%xmm3, %xmm1
	movdqa	%xmm2, %xmm3
	punpckldq	%xmm1, %xmm3
	punpckhdq	%xmm1, %xmm2
	movups	%xmm3, (%rdi,%rax)
	movups	%xmm2, 16(%rdi,%rax)
	addq	$32, %rax
	cmpq	$1024, %rax
	jne	.L307
	ret
.L309:
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	xorl	%edx, %edx
	pushq	%rbx
	.cfi_def_cfa_offset 24
	.cfi_offset 3, -24
	.p2align 4,,10
	.p2align 3
.L306:
	movl	(%r8,%rdx), %ebx
	movl	(%rsi,%rdx), %ebp
	movl	4(%r8,%rdx), %ecx
	movq	%rbx, %rax
	imulq	%rbp, %rax
	movq	%rcx, %r11
	imulq	%rbp, %rcx
	movl	%eax, %r10d
	shrq	$31, %rax
	andl	$2147483647, %r10d
	addl	%eax, %r10d
	movl	4(%rsi,%rdx), %eax
	imulq	%rax, %r11
	imulq	%rbx, %rax
	movl	%r11d, %r9d
	shrq	$31, %r11
	andl	$2147483647, %r9d
	addl	%r11d, %r9d
	movl	%ecx, %r11d
	shrq	$31, %rcx
	andl	$2147483647, %r11d
	addl	%r11d, %ecx
	movl	%eax, %r11d
	shrq	$31, %rax
	andl	$2147483647, %r11d
	addl	%eax, %r11d
	movl	%r9d, %eax
	shrl	$31, %r9d
	andl	$2147483647, %eax
	addl	%r9d, %eax
	movl	%r10d, %r9d
	shrl	$31, %r10d
	andl	$2147483647, %r9d
	notl	%eax
	addl	%r10d, %r9d
	andl	$2147483647, %eax
	addl	%r9d, %eax
	movl	%eax, %r9d
	shrl	$31, %eax
	andl	$2147483647, %r9d
	addl	%r9d, %eax
	movl	%ecx, %r9d
	shrl	$31, %ecx
	movl	%eax, (%rdi,%rdx)
	movl	%r11d, %eax
	andl	$2147483647, %r9d
	andl	$2147483647, %eax
	addl	%r9d, %eax
	addl	%ecx, %eax
	shrl	$31, %r11d
	addl	%r11d, %eax
	movl	%eax, %ecx
	shrl	$31, %eax
	andl	$2147483647, %ecx
	addl	%ecx, %eax
	movl	%eax, 4(%rdi,%rdx)
	addq	$8, %rdx
	cmpq	$1024, %rdx
	jne	.L306
	popq	%rbx
	.cfi_def_cfa_offset 16
	popq	%rbp
	.cfi_def_cfa_offset 8
	ret
	.cfi_endproc
.LFE51:
	.size	mrsn_mulc_256, .-mrsn_mulc_256
	.globl	omegas128
	.data
	.align 32
	.type	omegas128, @object
	.size	omegas128, 512
omegas128:
	.long	1
	.long	0
	.long	0
	.long	1
	.long	32768
	.long	32768
	.long	2147450879
	.long	32768
	.long	1556715293
	.long	978592373
	.long	1168891274
	.long	1556715293
	.long	978592373
	.long	1556715293
	.long	590768354
	.long	978592373
	.long	1241207368
	.long	1179735656
	.long	967747991
	.long	1241207368
	.long	2112881577
	.long	1415090252
	.long	732393395
	.long	2112881577
	.long	1415090252
	.long	2112881577
	.long	34602070
	.long	1415090252
	.long	1179735656
	.long	1241207368
	.long	906276279
	.long	1179735656
	.long	1641940819
	.long	26164677
	.long	2121318970
	.long	1641940819
	.long	1690787918
	.long	579625837
	.long	1567857810
	.long	1690787918
	.long	1133522282
	.long	280947147
	.long	1866536500
	.long	1133522282
	.long	567259857
	.long	194696271
	.long	1952787376
	.long	567259857
	.long	194696271
	.long	567259857
	.long	1580223790
	.long	194696271
	.long	280947147
	.long	1133522282
	.long	1013961365
	.long	280947147
	.long	579625837
	.long	1690787918
	.long	456695729
	.long	579625837
	.long	26164677
	.long	1641940819
	.long	505542828
	.long	26164677
	.long	206059115
	.long	1935040570
	.long	212443077
	.long	206059115
	.long	1796741361
	.long	1263730590
	.long	883753057
	.long	1796741361
	.long	408478793
	.long	262191051
	.long	1885292596
	.long	408478793
	.long	373229752
	.long	1309288441
	.long	838195206
	.long	373229752
	.long	228509164
	.long	14530030
	.long	2132953617
	.long	228509164
	.long	134155457
	.long	1038945916
	.long	1108537731
	.long	134155457
	.long	2079025011
	.long	2137679949
	.long	9803698
	.long	2079025011
	.long	2140339328
	.long	1742797653
	.long	404685994
	.long	2140339328
	.long	1742797653
	.long	2140339328
	.long	7144319
	.long	1742797653
	.long	2137679949
	.long	2079025011
	.long	68458636
	.long	2137679949
	.long	1038945916
	.long	134155457
	.long	2013328190
	.long	1038945916
	.long	14530030
	.long	228509164
	.long	1918974483
	.long	14530030
	.long	1309288441
	.long	373229752
	.long	1774253895
	.long	1309288441
	.long	262191051
	.long	408478793
	.long	1739004854
	.long	262191051
	.long	1263730590
	.long	1796741361
	.long	350742286
	.long	1263730590
	.long	1935040570
	.long	206059115
	.long	1941424532
	.long	1935040570
	.globl	omegas512
	.align 32
	.type	omegas512, @object
	.size	omegas512, 504
omegas512:
	.long	430821412
	.long	1152650470
	.long	236104903
	.long	1577470940
	.long	485600145
	.long	224958826
	.long	206059115
	.long	1935040570
	.long	660017901
	.long	1340846354
	.long	1896945393
	.long	2098580229
	.long	735494074
	.long	1494204761
	.long	1641940819
	.long	26164677
	.long	1777644782
	.long	1383853684
	.long	1093071961
	.long	648593218
	.long	1920912571
	.long	914097328
	.long	1742797653
	.long	2140339328
	.long	1371669334
	.long	2103108137
	.long	1260750973
	.long	1362440376
	.long	1336950523
	.long	839591040
	.long	1241207368
	.long	1179735656
	.long	1202912605
	.long	1980032781
	.long	1921627098
	.long	1668363411
	.long	165851886
	.long	1674906685
	.long	228509164
	.long	14530030
	.long	1792244284
	.long	36557796
	.long	1014093253
	.long	2137011181
	.long	252929270
	.long	1353673049
	.long	194696271
	.long	567259857
	.long	853979252
	.long	1113159341
	.long	1563928157
	.long	849605071
	.long	472916039
	.long	952794586
	.long	1309288441
	.long	373229752
	.long	528066207
	.long	951582730
	.long	141956360
	.long	1977033713
	.long	605970061
	.long	1530121874
	.long	1556715293
	.long	978592373
	.long	378535762
	.long	1207781610
	.long	477953613
	.long	2022380190
	.long	1398069297
	.long	1397384897
	.long	408478793
	.long	262191051
	.long	1326503162
	.long	1362518885
	.long	1276547035
	.long	1514613395
	.long	328267072
	.long	484667533
	.long	1133522282
	.long	280947147
	.long	1730666434
	.long	2110668387
	.long	81378258
	.long	1357626641
	.long	655387905
	.long	1395419301
	.long	1038945916
	.long	134155457
	.long	1603661239
	.long	1398285837
	.long	1949783546
	.long	1067683608
	.long	1916124599
	.long	187158958
	.long	1415090252
	.long	2112881577
	.long	834535867
	.long	222141861
	.long	1553669210
	.long	736262640
	.long	1925205788
	.long	278287463
	.long	2079025011
	.long	2137679949
	.long	1756768506
	.long	1678097410
	.long	1506666447
	.long	445356670
	.long	1370602608
	.long	1664948088
	.long	579625837
	.long	1690787918
	.long	1082787046
	.long	2133873350
	.long	1895558694
	.long	636875771
	.long	1644164930
	.long	820860779
	.long	1263730590
	.long	1796741361
	.long	578660954
	.long	152276873
	.long	2085743640
	.long	812986380
	.long	1854234209
	.long	1637799161
	.section	.rodata.cst16,"aM",@progbits,16
	.align 16
.LC0:
	.long	2147483647
	.long	2147483647
	.long	2147483647
	.long	2147483647
	.set	.LC1,.LC0
	.align 16
.LC3:
	.long	2147450880
	.long	2147450880
	.long	2147450880
	.long	2147450880
	.set	.LC4,.LC3
	.ident	"GCC: (Ubuntu 13.3.0-6ubuntu2~24.04.1) 13.3.0"
	.section	.note.GNU-stack,"",@progbits
	.section	.note.gnu.property,"a"
	.align 8
	.long	1f - 0f
	.long	4f - 1f
	.long	5
0:
	.string	"GNU"
1:
	.align 8
	.long	0xc0000002
	.long	3f - 2f
2:
	.long	0x3
3:
	.align 8
4:
