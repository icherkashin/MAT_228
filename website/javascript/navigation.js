var menuNodeArray = document.getElementsByClassName("menuNode")

      function display(event) {
      
	event.preventDefault() 
      
      if (this.nextSibling.nextSibling.style.visibility == "visible") 
	{ 
	  this.nextSibling.nextSibling.style.transition = "max-height 0.25s ease, opacity 0.25s ease, visibility 0.25s ease"	  
	  this.nextSibling.nextSibling.style.opacity = "0"
	  this.nextSibling.nextSibling.style.maxHeight = "0"
	  this.nextSibling.nextSibling.style.visibility = "hidden"
	} 
	else 
	{ 	
	  this.nextSibling.nextSibling.style.transition = "max-height 0.5s ease, opacity 1.25s ease, visibility 1.25s ease"	     
	  this.nextSibling.nextSibling.style.maxHeight = "253px"	
	  this.nextSibling.nextSibling.style.opacity = "1"
	  this.nextSibling.nextSibling.style.visibility = "visible"
	}
      }
           
      for (var i=0; i<menuNodeArray.length; i++) {
	  menuNodeArray[i].addEventListener("click", display);
	}	 
