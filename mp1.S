
.data					# section declaration

        # Useful offset constants for accessing members of a 
        # struct mp1_blink_struct structure
        LOCATION   = 0    # w
        ON_CHAR    = 2    # b
        OFF_CHAR   = 3    # b
        ON_LENGTH  = 4    # w
        OFF_LENGTH = 6    # w
        COUNTDOWN  = 8    # w
        STATUS     = 10   # w
        NEXT       = 12   # l

        STRUCT_SIZE = 16

# Pointer to head of list (initialized to NULL)
mp1_list_head:
        .long   0

.text					# section declaration

# Export the function symbol names

.global mp1_rtc_tasklet
.global mp1_ioctl

# void mp1_poke(void);
#
# Interface: Register-based arguments (not C-style)
#    Inputs: %cl  - The byte you wish to write
#            %eax - Offset from the start of video memory that you wish
#                   to write to
#   Outputs: Text-mode video screen is written to at location %eax with
#            the byte in %cl
# Registers: Clobbers EDX
mp1_poke:
        
	movl    vmem_base_addr(,1),%edx
	movb    %cl,(%edx,%eax,1)
	ret

mp1_rtc_tasklet:
	ret

mp1_ioctl:
        pushl %ebp
        movl %esp,%ebp
        movl 12(%ebp), %eax # fetch arg cmd to %eax
        # valid cmd test
        cmpl $0,%eax  # cmd>=0
        jl ioctl_invalid_cmd
        cmpl $3,%eax  # cmd<=3
        jg ioctl_invalid_cmd
        jmp *jump_table(,%eax,4) 
ioctl_invalid_cmd:
        movl $-1,%eax
        leave
        ret

mp1_ioctl_add:
        # call mp1_malloc: allocate mem for copying the struct to kernal
        pushl $16            # arg size
        call mp1_malloc     
        addl $4,%esp         # pop args 
        # %eax = to pointer
        cmpl $0,%eax
        je fail_mem_allowcate 
        pushl %eax           # save "to pointer" to local:-4(%ebp) 
        movl 8(%ebp), %ecx   # fetch arg "from pointer" to %ecx
        # check if "from pointer" is NULL
        cmpl $0, %ecx
        je fail_mem_copy
        # call mp1_copy_from_user: copy struct from user mem to kernal mem
        pushl $16            # arg size
        pushl %ecx           # arg "from pointer"
        pushl %eax           # arg "to pointer"
        call mp1_copy_from_user
        addl $12,%esp        # pop args 
        # %eax = number of bytes fail to copy
        cmpl $0,%eax
        jne fail_mem_copy
        # check the location is valid
        movl -4(%ebp), %eax  # load "to pointer" to %eax
        movw LOCATION(%eax), %cx  # load "location" to %cx
        cmpw $1999, %cx     # invalid "location" >1999=80*25-1
        jg fail_invalid_location
        cmpw $0, %cx        # invalid "location" <0
        jl fail_invalid_location
        pushw %cx            # save "location" to local:-6(%ebp) 
        # set the countdown field to be equal to the on length field
        movw ON_LENGTH(%eax),%cx  # fetch "ON_LENGTH" to cx
        movw %cx, COUNTDOWN(%eax)
        # set the status field to 1
        movw $1, STATUS(%eax)
        # add this item to the linked list
        movl mp1_list_head, %ecx # fetch "head's pointer" to ecx
        movl %ecx, NEXT(%eax)
        movl %eax, mp1_list_head
        # call mp1 poke
        movl %eax, %edx     # fetch "to pointer" to %edx
        xorl %eax, %eax
        movw -6(%ebp), %ax  # load "location" to ax
        addl %eax,%eax # doulbe eax
        movb ON_CHAR(%edx),%cl
        call mp1_poke
        # return 0 
        xorl %eax, %eax       
        leave
        ret

fail_mem_allowcate:
        # return -1
        movl $-1,%eax
        leave
        ret
fail_invalid_location:
fail_mem_copy:
        # free the mem by "to pointer"
        movl -4(%ebp), %eax  # load "to pointer" to %eax
        # call mp1_free: allocate mem for copying the struct to kernal
        pushl %eax            # argument to pointer
        call mp1_free     
        addl $4,%esp         # pop args 
        # return -1
        movl $-1,%eax
        leave
        ret 


        
mp1_ioctl_remove:
	ret
        
mp1_ioctl_find:
	ret

mp1_ioctl_sync:
	ret
/* The jump table for mp1_ioctl */
jump_table:
.long mp1_ioctl_add, mp1_ioctl_remove, mp1_ioctl_find, mp1_ioctl_sync

.end
