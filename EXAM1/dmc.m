% DMC
% klasa implementujaca algorytm DMC
classdef DMC <handle
   properties
    % dlugosc symulacji   
    simulation_length
    % odpowiedz skokowa
    s

    % horyzont dynamiki
    D
    % horyzont predykcji
    N
    % horyzont sterowania
    Nu
    % wspoczynniki wagowe
    lambda
    psi
    
    % liczba wyjsc
    ny
    % liczba wejsc
    nu
    
    
    K
    K1
   
    M
    Mp
    dUp
    
    dU

    Y_zad=0;
    E;
    
    e
    u
    y

   end
   methods
     
     function obj=DMC(D,N,Nu,lambda,psi,ny,nu,response)
        obj.s=response;
        obj.D=D;
        obj.N=N;
        obj.Nu=Nu;
    
        % z hardcodowane
        obj.ny=ny;
        obj.nu=nu;
        
        % wypelnianie macierzy
        obj.M=obj.fillM()
        obj.Mp=obj.fillMp()
        
        obj.psi=eye(obj.N*obj.ny,obj.N*obj.ny).*repmat(psi,obj.N,1);
        obj.lambda=eye(obj.Nu*obj.nu,obj.Nu*obj.nu).*repmat(lambda,obj.Nu,1);
    
        % deklaracja macierzy K
        obj.K=inv(obj.M'*obj.psi*obj.M+obj.lambda)*obj.M'*obj.psi;
        obj.K1=obj.K(1:obj.nu,:);
           
        obj.dUp=zeros((obj.D-1)*obj.nu,1);
     end
    
    
    function M=fillM(obj)
    M = zeros(obj.N*obj.ny,obj.Nu*obj.nu);
    for j=1:obj.Nu 
        for i=1:obj.N
            if(i>=j)
                M(1+(i-1)*obj.ny:i*obj.ny,1+(j-1)*obj.nu:j*obj.nu)=obj.s(i-j+1,:,:);
            end
        end
    end
    end

    % uzupełnianie macierzy Mp
    function Mp=fillMp(obj)
        Mp=zeros(obj.N*obj.ny,(obj.D-1)*obj.nu);
        for i=1:obj.N
            for j=1:obj.D-1
                if i+j<=obj.D
                    Mp(1+(i-1)*obj.ny:i*obj.ny,1+(j-1)*obj.nu:j*obj.nu)=obj.s(i+j,:,:)-obj.s(j,:,:);
                else
                    Mp(1+(i-1)*obj.ny:i*obj.ny,1+(j-1)*obj.nu:j*obj.nu)=obj.s(obj.D,:,:)-obj.s(j,:,:);
                end
            end
        end
    end
    
    function setup_simulation(obj,Y_zad)
        obj.Y_zad=Y_zad;
        obj.simulation_length=length(Y_zad(:,1));

        obj.e=zeros(obj.simulation_length,obj.ny);
        obj.u=zeros(obj.simulation_length,obj.nu);
        obj.y=zeros(obj.simulation_length,obj.ny);
    
    end

    function E=zadanie5(obj,variable)
        variable
        psi=variable(1:3)';
        lambda=variable(4:7)';
        obj.psi=eye(obj.N*obj.ny,obj.N*obj.ny).*repmat(psi,obj.N,1);
        obj.lambda=eye(obj.Nu*obj.nu,obj.Nu*obj.nu).*repmat(lambda,obj.Nu,1);
    
        % deklaracja macierzy K
        obj.K=inv(obj.M'*obj.psi*obj.M+obj.lambda)*obj.M'*obj.psi;
        obj.K1=obj.K(1:obj.nu,:);
        obj.simulation_process();
        E=sum(obj.E)
    end
    
    function simulation_process(obj)
       for k=5:obj.simulation_length

           [obj.y(k,1),obj.y(k,2),obj.y(k,3)]=symulacja_obiektu10(  obj.u(k-1,1),obj.u(k-2,1),obj.u(k-3,1),obj.u(k-4,1),...
                                                                    obj.u(k-1,2),obj.u(k-2,2),obj.u(k-3,2),obj.u(k-4,2),...
                                                                    obj.u(k-1,3),obj.u(k-2,3),obj.u(k-3,3),obj.u(k-4,3),...
                                                                    obj.u(k-1,4),obj.u(k-2,4),obj.u(k-3,4),obj.u(k-4,4),...
                                                                    obj.y(k-1,1),obj.y(k-2,1),obj.y(k-3,1),obj.y(k-4,1),...
                                                                    obj.y(k-1,2),obj.y(k-2,2),obj.y(k-3,2),obj.y(k-4,2),...
                                                                    obj.y(k-1,3),obj.y(k-2,3),obj.y(k-3,3),obj.y(k-4,3));
            % wyliczenie uchybu
            obj.e(k,:)=obj.Y_zad(k,:)-obj.y(k,:);
           
            % wyliczanie zmiany sterowania
            obj.dU=obj.K1*(repmat(obj.e(k,:)',obj.N,1)-obj.Mp*obj.dUp);

            % wyznaczanie sterowania
            obj.u(k,:)=obj.u(k-1,:)+obj.dU';

            % wyznaczanie wektora poprzednich delt sterowania
            obj.dUp=circshift(obj.dUp,obj.nu);
            obj.dUp(1:obj.nu)=(obj.u(k,:)-obj.u(k-1,:))';
        end
        obj.E=sum(obj.e.^2);
    end    
    
    function simulation_classic(obj)
       for k=5:obj.simulation_length

           [obj.y(k,1),obj.y(k,2),obj.y(k,3)]=symulacja_obiektu10(  obj.u(k-1,1),obj.u(k-2,1),obj.u(k-3,1),obj.u(k-4,1),...
                                                                    obj.u(k-1,2),obj.u(k-2,2),obj.u(k-3,2),obj.u(k-4,2),...
                                                                    obj.u(k-1,3),obj.u(k-2,3),obj.u(k-3,3),obj.u(k-4,3),...
                                                                    obj.u(k-1,4),obj.u(k-2,4),obj.u(k-3,4),obj.u(k-4,4),...
                                                                    obj.y(k-1,1),obj.y(k-2,1),obj.y(k-3,1),obj.y(k-4,1),...
                                                                    obj.y(k-1,2),obj.y(k-2,2),obj.y(k-3,2),obj.y(k-4,2),...
                                                                    obj.y(k-1,3),obj.y(k-2,3),obj.y(k-3,3),obj.y(k-4,3));
            % wyliczenie uchybu
            obj.e(k,:)=obj.Y_zad(k,:)-obj.y(k,:);
           
            % wyliczanie zmiany sterowania
            obj.dU=obj.K*(repmat(obj.e(k,:)',obj.N,1));

            % wyznaczanie sterowania
            obj.u(k,:)=obj.u(k-1,:)+obj.dU(1:4)';%??????????????????

            % wyznaczanie wektora poprzednich delt sterowania
            obj.dUp=circshift(obj.dUp,obj.nu);
            obj.dUp(1:obj.nu)=(obj.u(k,:)-obj.u(k-1,:))';
        end
        obj.E=sum(obj.e.^2);
    end    
end
end