#ifndef _TEXTFORRENDER_H_
#define _TEXTFORRENDER_H_

#include <string>
#include <set>
#include <boost/shared_ptr.hpp>
#include <QFont>

namespace QGLVEXT{

  /**
   * @class TextWithPosition string type text with position information.
   * 
   */
  class TextWithPosition{
	
  public:
	TextWithPosition(const std::string text="",int x=1,int y=1,const QFont font = QFont("Times", 30, QFont::Bold)){
	  setContent(text);
	  setX(x);
	  setY(y);
	  setFont(font);
	}
	TextWithPosition(const TextWithPosition &other){
	  *this = other;
	}

	TextWithPosition & operator = (const TextWithPosition &other){
	  _content = other.getContent();
	  _x = other.getX();
	  _y = other.getY();
	  _font = other.getFont();
	  return (*this);
	}

	void setContent(const std::string content){
	  _content = content;
	}
	void setFont(const QFont &font){
	  _font = font;
	}
	void setX(const int x){
	  _x = x;
	}
	void setY(const int y){
	  _y = y;
	}

	const std::string &getContent()const{
	  return _content;
	}
	const QFont &getFont()const{
	  return _font;
	}
	int getX()const{
	  return _x;
	}
	int getY()const{
	  return _y;
	}
	
  protected:
	int _x;
	int _y;
	std::string _content;
	QFont _font;
  };

  typedef boost::shared_ptr< TextWithPosition > pTextWithPosition; 

  /**
   * @class TextWithPositionCompare compare two TextWithPosition elements, used
   * for define TextPosVec bellow.
   * 
   */
  class TextWithPositionCompare{
  public:
	bool operator() (const TextWithPosition &lhs, const TextWithPosition& rhs) const{
	  const int x1 = lhs.getX();
	  const int y1 = lhs.getY();
	  const int x2 = rhs.getX();
	  const int y2 = rhs.getY();
	  return (x1<x2)||(x1==x2&&y1<y2);
	}
  };
  
  typedef std::set<TextWithPosition,TextWithPositionCompare> TextSet;

  class TextForRender: public TextSet{
	
  public:
    void update(const TextWithPosition &text){
	  TextSet::erase(text);
	  TextSet::insert(text);
	}
	void update(const std::string content, int x, int y){
	  const TextWithPosition text(content,x,y);
	  TextSet::erase(text);
	  TextSet::insert(text);
	}
  };
  typedef boost::shared_ptr<TextForRender > pTextForRender;
  
}//end of namespace

#endif /*_TEXTFORRENDER_H_*/
