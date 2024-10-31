#ifndef STACK_H
#define STACK_H

#include <assert.h>
#include <stdio.h>

//-----------------------------------------------------------------------
//
// 
//-----------------------------------------------------------------------
template<class T>
class STACK 
{
	protected:
		T*  m_buffer     ;
		int m_size       ;
		int m_size_buffer;

	public:
		//-----------------------------------------------------------------------
		// constructor
		//-----------------------------------------------------------------------
		STACK() 
		{
			m_size=0;
			m_size_buffer = 100;
			m_buffer = new T[m_size_buffer];
		}

		STACK( const STACK& st) 
		{
			m_size=st.m_size;
			m_size_buffer = st.m_size_buffer;
			m_buffer = new T[m_size_buffer];
			memcpy(m_buffer,st.m_buffer,sizeof(T)*m_size);
		}

		~STACK()
		{
			delete[] m_buffer;
		}

		void operator=( const STACK& St )
		{
			delete[] m_buffer;
			m_size        = St.m_size;
			m_size_buffer = St.m_size_buffer;
			m_buffer      = new T[m_size_buffer];

			for(int n=0;n<m_size;n++)
				m_buffer[n] = St.m_buffer[n];
		}

		//-----------------------------------------------------------------------
		// accessors
		//-----------------------------------------------------------------------
		const T& operator[](int i) const{ return m_buffer[i]; }
		const T& operator()(int i) const{ return m_buffer[i]; }
		      T& operator[](int i)      { return m_buffer[i]; }
		      T& operator()(int i)      { return m_buffer[i]; }

		//-----------------------------------------------------------------------
		// conversion
		//-----------------------------------------------------------------------
		operator const T*() const{ return (const T*)m_buffer;}
		operator       T*() const{ return (      T*)m_buffer;}

		//-----------------------------------------------------------------------
		// i = 0,1,...,size-1
		//-----------------------------------------------------------------------
		int size() const { return m_size;}
		bool is_empty() const { return m_size==0; }

		//-----------------------------------------------------------------------
		// print out
		//-----------------------------------------------------------------------
		void print(const char* name) const
		{
			printf("-----------------STACK : %s---------------------\n",name);

			for(int n=0;n<m_size;n++)
			{
				if(sizeof(T)==sizeof(double)) printf("%d th : %f\n",n, m_buffer[n]);
				if(sizeof(T)==sizeof(int   )) printf("%d th : %d\n",n, m_buffer[n]);
				if(sizeof(T)==sizeof(char  )) printf("%d th : %c\n",n, m_buffer[n]);
			}
		}

		//-----------------------------------------------------------------------
		// Push : add    the element at the end
		// Pop  : remove the element at the end
		//-----------------------------------------------------------------------
		void push( const T& element )
		{
			if(m_size_buffer==m_size)
			{
				m_size_buffer = m_size + (int)(m_size*.5);
				
				T* new_buffer = new T[m_size_buffer];
				if(new_buffer==0){
					printf("memory allocation fails in STACK\n");
					assert(false);
				}
				for(int n=0;n<m_size;n++) new_buffer[n] = m_buffer[n];

				delete[] m_buffer; 
				m_buffer = new_buffer;
			}

			m_buffer[m_size] = element;
			m_size++;
		}

		T pop()
		{
			if(m_size==0) assert(false);

			T out = m_buffer[m_size-1];
			m_size--;
			return out;
		}

		void make_it_empty()
		{
			m_size = 0;
		}
};

#endif
